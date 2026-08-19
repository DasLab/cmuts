/* pipeline.c -- the loader, the worker pool, and the completion consumer.
 *
 * The input is read once, sequentially, by a single loader. Each read is copied
 * into a pooled carrier and queued; workers take them in batches, accumulate
 * into a private shadow, and merge into the shared per-reference accumulator
 * only when they move to another reference. A reference is finished when its
 * last handle is dropped, at which point it goes to the consumer to be written
 * and recycled.
 *
 * Every buffer is drawn from a fixed pool, so memory is bounded by the pool
 * sizes and the longest reference, not by the size of the file.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "pipeline.h"

#include <stdatomic.h>

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bamstream.h"
#include "h5writer.h"
#include "itempool.h"
#include "progress.h"
#include "queue.h"
#include "refctx.h"
#include "refrow.h"
#include "refseq.h"
#include "tally.h"

/* Completed references drained per pop. The consumer is never the bottleneck, so
 * this only saves a lock acquisition per reference. */
#define COMPLETION_BATCH 16

#define DEFAULT_WORKERS        1
#define DEFAULT_DECODE_THREADS 0
#define DEFAULT_QUEUE_CAPACITY 4096
#define DEFAULT_BATCH          64
#define DEFAULT_LIVE_REFS      64

/* Output chunks offered to the workers and not yet written back. The bound is what lets
 * the workers hand filtered chunks over without ever blocking: the filtered queue holds
 * this many, and the consumer offers no more until some come back. */
#define CHUNKS_PER_WORKER 2

typedef struct {
    cm_bam_stream *bam;
    refseq_source *refs;
    queue         *work;       /* loader -> workers */
    queue         *completed;  /* last handle dropped -> consumer */
    queue         *filtered;   /* chunks the workers filtered -> consumer */
    size_t         chunk_cap;  /* chunks out with the workers at once */
    itempool      *items;
    ctxpool       *contexts;
    h5writer      *out;
    refrow        *rows;       /* writes a finished accumulator */
    progress      *bar;
    tally_tables   tally_tables;
    filter_config  filter_config;
    size_t         batch;
    size_t         ref_cap;    /* longest reference, sizing every accumulator */
    bool           pairwise;
    bool           wanted[OUT_N_FIELDS];   /* the optional fields this run writes */
} pipeline;

/* The pipeline is filled in while the run is assembled and only read once it starts,
 * so every stage below takes it const. */

/* The first failure a worker reported. One is enough to stop the run.
 *
 * A struct so that the atomic is reachable only through the two calls below. Assigning
 * straight to it would keep the last failure offered rather than the first. */
typedef struct {
    _Atomic int status;   /* a phmm_status, PHMM_OK until a worker stops */
} failure_flag;

/* Records a failure, keeping the first offered and discarding the rest. */
static void failure_record(failure_flag *f, phmm_status status)
{
    int unfailed = PHMM_OK;

    atomic_compare_exchange_strong(&f->status, &unfailed, (int)status);
}

static phmm_status failure_seen(const failure_flag *f)
{
    return (phmm_status)atomic_load(&f->status);
}

/* What each way a marginalization can fail is reported as. PHMM_NO_PATH is listed for the
 * switch and not reached: the worker counts that read as rejected and reads on. */
static const char *failure_text(phmm_status status)
{
    switch (status) {
        case PHMM_NO_MEMORY:
            return "out of memory marginalizing a read";
        case PHMM_NO_PATH:
            return "no alignment of a read has any probability under these rates";
        case PHMM_UNSOUND:
            return "a marginalization did not hold together";
        case PHMM_OK:
            break;
    }

    return "a read could not be marginalized";
}

static void pipeline_finish_reference(const pipeline *p, refctx *ctx)
{
    void *handle = ctx;

    queue_push_all(p->completed, &handle, 1);
}

/* ------------------------------------------------------------------------ */
/* Workers                                                                   */
/* ------------------------------------------------------------------------ */

typedef struct {
    const pipeline *pipe;
    failure_flag  *failure;  /* shared; the one thing a worker writes */
    refctx        *held;    /* reference the shadow holds a handle for, or NULL */
    accum          shadow;
    pairs          shadow_pairs;   /* cells only under --pairwise */
    tally_scratch *scratch; /* buffers for the processing step */
    void         **slots;   /* a batch of workitems, allocated before the
                               thread starts */
    pthread_t      thread;
} worker;

/* Merges the shadow into its reference and releases the handle protecting it. The
 * shadow is left entirely zero. */
static void worker_flush_shadow(worker *w)
{
    refctx *ctx = w->held;

    if (!ctx) {
        return;
    }

    refctx_merge(ctx, &w->shadow, w->shadow_pairs.cells ? &w->shadow_pairs : NULL);
    accum_zero(&w->shadow, ctx->len);

    if (w->shadow_pairs.cells) {
        pairs_zero(&w->shadow_pairs, ctx->len);
    }

    w->held = NULL;

    if (refctx_release(ctx, 1)) {
        pipeline_finish_reference(w->pipe, ctx);
    }
}

static void worker_switch_shadow(worker *w, refctx *ctx)
{
    if (w->held == ctx) {
        return;
    }

    worker_flush_shadow(w);
    refctx_acquire(ctx, 1);
    w->held = ctx;
}

/* Returns one slot of a batch as the workitem it is. The queue and the pool carry void *,
 * so the conversion back is written here once. */
static const workitem *item_at(void *const *slots, size_t i)
{
    return slots[i];
}

/* Counts every read of a run into the worker's own shadow. */
static void worker_count_run(worker *w, void **slots, size_t n,
                             const cm_fasta_record *ref)
{
    for (size_t i = 0; i < n; i++) {
        const workitem *item = item_at(slots, i);
        cm_bam_record   read;
        phmm_status     status;

        cm_bam_record_view(item->rec, &read);
        status = tally(&read, ref, &w->pipe->tally_tables, w->scratch,
                       &w->shadow,
                       w->shadow_pairs.cells ? &w->shadow_pairs : NULL);

        /* The tally has counted a read it could not score as rejected, so the run goes
         * on. Every other failure would meet every read after it. */
        if (status == PHMM_NO_PATH) {
            continue;
        }

        if (status != PHMM_OK) {
            failure_record(w->failure, status);
            return;
        }
    }
}

/* Processes a run of reads belonging to one reference. */
static void worker_process_run(worker *w, void **slots, size_t n)
{
    refctx         *ctx   = item_at(slots, 0)->ctx;
    cm_fasta_record ref;

    worker_switch_shadow(w, ctx);
    refctx_sequence(ctx, &ref);

    /* Once a worker has failed the run is released without being counted, since the
     * output will never be written. */
    if (failure_seen(w->failure) == PHMM_OK) {
        worker_count_run(w, slots, n, &ref);
    }

    /* The whole run's carriers go back at once. Nothing reads them afterwards:
     * worker_process_batch scans only at or beyond the current head. */
    itempool_give_many(w->pipe->items, slots, n);

    /* The shadow took a handle on this reference before the run was counted, so this is
     * never the last one and the reference cannot finish here. */
    refctx_release(ctx, (int)n);
}

/* Filters an output chunk on behalf of the writing thread, which is how the compression
 * shares the threads the run was given. The carrier was allocated for this one trip. */
static void worker_filter_chunk(const pipeline *p, workitem *item)
{
    void *chunk = item->chunk;

    h5chunk_filter(item->chunk);
    queue_push_all(p->filtered, &chunk, 1);
    free(item);
}

static void worker_process_batch(worker *w, void **slots, size_t n)
{
    for (size_t i = 0; i < n; ) {
        const workitem *head = slots[i];
        size_t          run  = 1;

        if (head->chunk) {
            worker_filter_chunk(w->pipe, slots[i]);
            i++;
            continue;
        }

        while (i + run < n && item_at(slots, i + run)->ctx == head->ctx) {
            run++;
        }

        worker_process_run(w, slots + i, run);
        i += run;
    }
}

static void *worker_main(void *arg)
{
    worker *w     = arg;
    size_t  batch = w->pipe->batch;

    for (;;) {
        size_t n = queue_try_pop(w->pipe->work, w->slots, batch);

        if (n == 0) {
            /* Flush before blocking: an idle worker holding a shadow pins its
             * reference indefinitely, and the loader then stalls waiting for a
             * context to come free. */
            worker_flush_shadow(w);

            n = queue_pop(w->pipe->work, w->slots, batch);
            if (n == 0) {
                break;
            }
        }

        worker_process_batch(w, w->slots, n);
    }

    worker_flush_shadow(w);
    return NULL;
}

/* ------------------------------------------------------------------------ */
/* Loader                                                                    */
/* ------------------------------------------------------------------------ */

/* The loader's state.
 *
 * The three pairs below are each two values meaningful only together: a reference
 * and its count of filtered reads, a batch and how much of it is filled, a
 * reservoir of carriers and how many remain. */
typedef struct {
    const pipeline     *pipe;
    const failure_flag *failure;  /* shared; the loader only reads it */

    refctx *reference;  /* the one being filled, or NULL */
    size_t  rejected;   /* its reads the filter turned away */

    void  **batch;      /* workitems bound for the workers */
    size_t  queued;

    void  **spare;      /* unspent workitem carriers drawn from the pool */
    size_t  held;

    int32_t owed;       /* the first reference not yet accounted for */
} loader;

static int loader_open(loader *l, const pipeline *p, const failure_flag *f)
{
    *l = (loader){ .pipe = p, .failure = f };

    l->batch = calloc(p->batch, sizeof *l->batch);
    l->spare = calloc(p->batch, sizeof *l->spare);

    if (l->batch && l->spare) {
        return 0;
    }

    free(l->batch);
    free(l->spare);
    return -1;
}

/* Queues the batch for the workers, taking a handle on the reference for each read so
 * that it cannot be finished while any of them is in transit. A read holds no handle
 * between being taken and reaching here. The loader's own handle covers that gap, taken
 * when the reference is opened. */
static void loader_dispatch(loader *l)
{
    if (l->queued == 0) {
        return;
    }

    refctx_acquire(l->reference, (int)l->queued);
    queue_push_all(l->pipe->work, l->batch, l->queued);
    l->queued = 0;
}

/* Releases the current reference: queued reads first, then its filtered count,
 * then the loader's own handle. The count must arrive before the handle is
 * dropped, or the reference could be written without it. */
static void loader_leave_reference(loader *l)
{
    refctx *ctx = l->reference;

    if (!ctx) {
        return;
    }

    loader_dispatch(l);

    if (l->rejected) {
        refctx_add_scalar(ctx, ACCUM_FILTERED, (double)l->rejected);
    }

    l->reference = NULL;
    l->rejected  = 0;

    if (refctx_release(ctx, 1)) {
        pipeline_finish_reference(l->pipe, ctx);
    }
}

static refctx *pipeline_open_reference(const pipeline *p, int32_t tid)
{
    const cm_fasta_record *seq = refseq_advance(p->refs, tid);
    if (!seq) {
        return NULL;
    }

    refctx *ctx = ctxpool_take(p->contexts);
    if (!ctx) {
        return NULL;
    }

    refctx_open(ctx, tid, cm_bam_stream_refname(p->bam, tid), seq);
    h5writer_expect(p->out, tid);
    return ctx;
}

/* Opens and closes a reference the reader passed over, so that its row is written unread.
 * A reference that received nothing is zero everywhere, which is what the accumulator a
 * pooled context carries already holds. */
static bool loader_emit_empty(loader *l, int32_t tid)
{
    refctx *ctx = pipeline_open_reference(l->pipe, tid);

    if (!ctx) {
        return false;
    }

    if (refctx_release(ctx, 1)) {
        pipeline_finish_reference(l->pipe, ctx);
    }

    return true;
}

/* Opens and closes every reference the reader has passed since the last one it stopped at,
 * which checks each against the headers and writes whatever its fields need. A reference
 * that needs nothing written costs the open and no more. */
static bool loader_account_through(loader *l, int32_t upto)
{
    while (l->owed < upto) {
        if (!loader_emit_empty(l, l->owed++)) {
            return false;
        }
    }

    return true;
}

/* Moves to the reference a read belongs to, unless it is already the current one.
 * Called before the filter is applied, so a rejected read still counts against its
 * reference: one whose reads were all filtered out is distinct from one that
 * received none. */
static bool loader_on_reference(loader *l, int32_t tid)
{
    if (l->reference && l->reference->tid == tid) {
        return true;
    }

    loader_leave_reference(l);

    if (!loader_account_through(l, tid)) {
        return false;
    }

    l->reference = pipeline_open_reference(l->pipe, tid);
    l->owed      = tid + 1;

    return l->reference != NULL;
}

/* Returns one carrier for the read just taken. Carriers are drawn a batch at a time and
 * spent one by one, so the pool's lock is taken once per batch and not once per
 * read. */
static workitem *loader_carrier(loader *l)
{
    if (l->held == 0) {
        l->held = itempool_take_many(l->pipe->items, l->spare, l->pipe->batch);
    }

    return l->held ? l->spare[--l->held] : NULL;
}

/* Copies the record just read and queues it for the workers. The copy is needed because
 * the reader overwrites its own record on the next advance. */
static int loader_take_read(loader *l)
{
    workitem *item = loader_carrier(l);

    if (!item) {
        return -1;
    }

    if (!bam_copy1(item->rec, cm_bam_stream_raw(l->pipe->bam))) {
        l->spare[l->held++] = item;
        return -1;
    }

    item->ctx             = l->reference;
    l->batch[l->queued++] = item;

    if (l->queued == l->pipe->batch) {
        loader_dispatch(l);
    }

    return 0;
}

/* Returns everything the loader still holds, whichever way the loop was left. */
static void loader_finish(loader *l)
{
    loader_leave_reference(l);

    if (l->held) {
        itempool_give_many(l->pipe->items, l->spare, l->held);
    }

    free(l->spare);
    free(l->batch);
}

/* Reads the file once, in order, dispatching what survives. */
static int loader_main(const pipeline *p, const failure_flag *f,
                       size_t *unmapped, char *error, size_t error_len)
{
    loader        l;
    cm_bam_record rec;
    int           status = CM_ITER_EOF;
    int           result = 0;

    if (loader_open(&l, p, f) < 0) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    while ((status = cm_bam_stream_next(p->bam, &rec)) == CM_ITER_OK) {
        progress_follow(p->bar);

        /* A worker has failed, so reading on gains nothing. Not an error here: the
         * worker reports it. */
        if (failure_seen(l.failure) != PHMM_OK) {
            break;
        }

        /* Unmapped reads align to no reference, so they are counted for the run as
         * a whole and go no further. */
        if (rec.flag & BAM_FUNMAP) {
            (*unmapped)++;
            continue;
        }

        if (!loader_on_reference(&l, rec.tid)) {
            const char *cause = refseq_error(p->refs);

            snprintf(error, error_len, "%s", cause ? cause : "no context available");
            result = -1;
            break;
        }

        if (!filter_accepts(&p->filter_config, &rec)) {
            l.rejected++;
            continue;
        }

        if (loader_take_read(&l) < 0) {
            snprintf(error, error_len, "unable to take a copy of an alignment");
            result = -1;
            break;
        }
    }

    /* Released before the references the reader never reached are opened. Holding it
     * while asking the pool for another deadlocks when the pool holds one context, which
     * is what --live-refs 1 sets. */
    loader_leave_reference(&l);

    /* The references the reader stopped short of received nothing, as did those it
     * passed over between reads. */
    if (result == 0 && status == CM_ITER_EOF
        && !loader_account_through(&l, cm_bam_stream_nref(p->bam))) {
        snprintf(error, error_len, "%s", refseq_error(p->refs));
        result = -1;
    }

    progress_follow(p->bar);
    loader_finish(&l);

    if (result == 0 && status == CM_ITER_ERROR) {
        snprintf(error, error_len, "%s", cm_bam_stream_error(p->bam));
        result = -1;
    }

    return result;
}

/* ------------------------------------------------------------------------ */
/* Consumer                                                                  */
/* ------------------------------------------------------------------------ */

typedef struct {
    const pipeline *pipe;
    int       status;     /* first write failure, if any */
    size_t    in_flight;  /* chunks with the workers, not yet written back */
    pthread_t thread;
} consumer;

static void consumer_fail(consumer *c)
{
    if (c->status == 0) {
        c->status = -1;
    }
}

/* Writes every chunk the workers have handed back so far. */
static void consumer_collect_chunks(consumer *c)
{
    void *slot;

    while (queue_try_pop(c->pipe->filtered, &slot, 1) == 1) {
        if (h5writer_write_chunk(c->pipe->out, slot) < 0) {
            consumer_fail(c);
        }

        c->in_flight--;
    }
}

/* Offers one chunk to the workers, or filters it here. The workers refuse it where their
 * queue is full or already closed, and the filtered queue has room promised only for
 * chunk_cap at a time; a run whose chunks outpace the workers is one whose consumer has
 * the spare time to filter its own. */
static void consumer_send_chunk(consumer *c, h5chunk *chunk)
{
    if (c->in_flight < c->pipe->chunk_cap) {
        workitem *item = malloc(sizeof *item);
        void     *slot = item;

        if (item) {
            *item = (workitem){ .chunk = chunk };

            if (queue_try_push(c->pipe->work, &slot, 1) == 1) {
                c->in_flight++;
                return;
            }

            free(item);
        }
    }

    h5chunk_filter(chunk);
    if (h5writer_write_chunk(c->pipe->out, chunk) < 0) {
        consumer_fail(c);
    }
}

static void consumer_send_settled(consumer *c)
{
    h5chunk *chunk;

    while ((chunk = h5writer_take_chunk(c->pipe->out))) {
        consumer_send_chunk(c, chunk);
    }
}

/* Gives one reference its row. A reference no read arrived on hands over no accumulator,
 * having accumulated nothing for a field to be derived from. */
static int consumer_row(consumer *c, refctx *ctx)
{
    const accum *acc = refctx_accumulated(ctx) ? &ctx->acc : NULL;

    return refrow_write(c->pipe->rows, ctx->tid, ctx->len, acc,
                        ctx->pr.cells ? &ctx->pr : NULL);
}

/* Writes one reference's row and reports it to the writer. The report is made whatever
 * the write did: a chunk settles only once every reference opened in it is accounted
 * for. */
static void consumer_write_reference(consumer *c, refctx *ctx)
{
    /* Keep draining after a failure, or the loader and workers block on a queue
     * nothing is emptying. */
    if (c->status == 0 && consumer_row(c, ctx) < 0) {
        consumer_fail(c);
    }

    if (h5writer_wrote(c->pipe->out, ctx->tid) < 0) {
        consumer_fail(c);
    }
}

/* The chunks still with the workers have all been filtered: the workers are joined
 * before the completed queue closes, and they drain their queue before exiting. */
static void consumer_finish_chunks(consumer *c)
{
    void *slot;

    while (c->in_flight > 0 && queue_pop(c->pipe->filtered, &slot, 1) == 1) {
        if (h5writer_write_chunk(c->pipe->out, slot) < 0) {
            consumer_fail(c);
        }

        c->in_flight--;
    }
}

static void *consumer_main(void *arg)
{
    consumer *c = arg;
    void     *slots[COMPLETION_BATCH];
    size_t    n;

    while ((n = queue_pop(c->pipe->completed, slots, COMPLETION_BATCH)) > 0) {
        for (size_t i = 0; i < n; i++) {
            refctx *ctx = slots[i];

            consumer_write_reference(c, ctx);
            ctxpool_give(c->pipe->contexts, ctx);
            consumer_collect_chunks(c);
            consumer_send_settled(c);
        }
    }

    consumer_finish_chunks(c);
    return NULL;
}

/* ------------------------------------------------------------------------ */
/* Assembly                                                                  */
/* ------------------------------------------------------------------------ */

pipeline_config pipeline_defaults(void)
{
    return (pipeline_config){
        .rate_config    = rate_defaults(),
        .verify         = REFSEQ_VERIFY_ALL,
        .workers        = DEFAULT_WORKERS,
        .decode_threads = DEFAULT_DECODE_THREADS,
        .queue_capacity = DEFAULT_QUEUE_CAPACITY,
        .batch          = DEFAULT_BATCH,
        .live_refs      = DEFAULT_LIVE_REFS,
        .filter_config  = filter_defaults(),
        .tally_config   = tally_defaults(),
    };
}

static void pipeline_teardown(pipeline *p)
{
    progress_finish(p->bar);
    refrow_destroy(p->rows);
    h5writer_close(p->out);
    ctxpool_destroy(p->contexts);
    itempool_destroy(p->items);
    queue_destroy(p->filtered);
    queue_destroy(p->completed);
    queue_destroy(p->work);
    refseq_close(p->refs);
    cm_bam_stream_close(p->bam);
}

static int pipeline_open_inputs(pipeline *p, const pipeline_config *cfg,
                                char *error, size_t error_len)
{
    const char *why;

    p->bam = cm_bam_stream_open(cfg->bam_paths, cfg->n_bams, cfg->fasta_path,
                                cfg->decode_threads);
    if (!p->bam) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    if (cm_bam_stream_error(p->bam)) {
        snprintf(error, error_len, "%s", cm_bam_stream_error(p->bam));
        return -1;
    }

    p->refs = refseq_open(cfg->fasta_path, p->bam, cfg->verify, &why);
    if (!p->refs) {
        snprintf(error, error_len, "%s: %s", cfg->fasta_path, why);
        return -1;
    }

    return 0;
}

static int pipeline_build_buffers(pipeline *p, const pipeline_config *cfg,
                                  char *error, size_t error_len)
{
    /* Enough carriers for a full queue, a batch per worker, and two batches for the
     * loader: one being filled and one in reserve. A short refill can leave it holding
     * part of each. */
    size_t carriers = cfg->queue_capacity + (cfg->workers + 2) * cfg->batch;

    p->batch         = cfg->batch;
    p->pairwise      = cfg->pairwise != 0;
    p->filter_config = cfg->filter_config;
    p->ref_cap       = (size_t)cm_bam_stream_max_reflen(p->bam);

    p->chunk_cap = cfg->workers * CHUNKS_PER_WORKER;
    p->work      = queue_create(cfg->queue_capacity);
    p->completed = queue_create(cfg->live_refs);
    p->filtered  = queue_create(p->chunk_cap);
    p->items     = itempool_create(carriers);
    p->contexts  = ctxpool_create(cfg->live_refs, p->ref_cap, cfg->pairwise != 0);

    if (!p->work || !p->completed || !p->filtered || !p->items || !p->contexts) {
        snprintf(error, error_len, "out of memory building the pipeline");
        return -1;
    }

    return 0;
}

static int pipeline_open_output(pipeline *p, const pipeline_config *cfg,
                                const char *program, bool may_replace, char *error,
                                size_t error_len)
{
    p->out = h5writer_create(cfg->output_path, program,
                             cm_bam_stream_nref(p->bam), p->ref_cap, may_replace,
                             p->wanted, true);
    if (!p->out) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    if (h5writer_error(p->out)) {
        snprintf(error, error_len, "%s: %s", cfg->output_path,
                 h5writer_error(p->out));
        return -1;
    }

    p->rows = refrow_create(p->out, cfg->rate_config, p->ref_cap, p->wanted);
    if (!p->rows) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    return 0;
}

/* Starts the workers, reporting through started how many are running so that the caller
 * joins exactly those. Every allocation a worker needs is made before any thread starts,
 * so running out of memory is reported instead of leaving a thread that exits at once. */
static int worker_start_all(worker *workers, size_t n, const pipeline *p,
                            failure_flag *f, size_t *started,
                            char *error, size_t error_len)
{
    *started = 0;

    for (size_t i = 0; i < n; i++) {
        workers[i].pipe    = p;
        workers[i].failure = f;
        workers[i].slots   = calloc(p->batch, sizeof *workers[i].slots);
        workers[i].scratch = tally_scratch_create();

        if (!workers[i].slots || !workers[i].scratch
            || accum_alloc(&workers[i].shadow, p->ref_cap) < 0
            || (p->pairwise && pairs_alloc(&workers[i].shadow_pairs, p->ref_cap) < 0)) {
            snprintf(error, error_len, "out of memory preparing worker %zu", i);
            return -1;
        }
    }

    for (size_t i = 0; i < n; i++) {
        if (pthread_create(&workers[i].thread, NULL, worker_main, &workers[i]) != 0) {
            snprintf(error, error_len, "unable to start worker %zu", i);
            return -1;
        }

        *started += 1;
    }

    return 0;
}

int pipeline_run(const pipeline_config *cfg, const char *program,
                 const out_manifest *writes, char *error, size_t error_len)
{
    pipeline     p           = { 0 };
    failure_flag failed      = { 0 };
    consumer     cons        = { 0 };
    bool         may_replace = false;
    size_t       unmapped    = 0;
    worker      *workers     = calloc(cfg->workers, sizeof *workers);
    size_t       started     = 0;
    int          status      = -1;

    if (!workers) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    /* The manifest says what a run of this program writes; the squares among it are
     * written only where they were asked for, and the pairwise coverage whenever any
     * statistic is. */
    out_selection(writes, p.wanted);
    p.wanted[OUT_PAIRWISE_CORRELATION] &= (cfg->pairwise & PAIRS_CORRELATION) != 0;
    p.wanted[OUT_PAIRWISE_CONDITIONAL] &= (cfg->pairwise & PAIRS_CONDITIONAL) != 0;
    p.wanted[OUT_PAIRWISE_COVERAGE]    &= cfg->pairwise != 0;

    if (h5writer_may_replace(cfg->output_path, cfg->overwrite, &may_replace,
                             error, error_len) < 0 ||
        pipeline_open_inputs(&p, cfg, error, error_len) < 0 ||
        pipeline_build_buffers(&p, cfg, error, error_len) < 0 ||
        pipeline_open_output(&p, cfg, program, may_replace, error, error_len) < 0) {
        goto done;
    }

    tally_tables_build(&p.tally_tables, &cfg->tally_config);

    /* Started last, so that nothing draws a bar and then fails before the first read. */
    p.bar = progress_start(p.bam);

    cons.pipe = &p;
    if (pthread_create(&cons.thread, NULL, consumer_main, &cons) != 0) {
        snprintf(error, error_len, "unable to start the consumer");
        goto done;
    }

    if (worker_start_all(workers, cfg->workers, &p, &failed, &started,
                         error, error_len) == 0) {
        status = loader_main(&p, &failed, &unmapped, error, error_len);
    }

    queue_close(p.work);
    for (size_t i = 0; i < started; i++) {
        pthread_join(workers[i].thread, NULL);
    }

    queue_close(p.completed);
    pthread_join(cons.thread, NULL);

    /* Checked once every worker has stopped, which catches a failure on the last batch
     * as well as one the loader saw in time to stop for. A loader that failed on its own
     * account has already reported something more specific. */
    if (status == 0 && failure_seen(&failed) != PHMM_OK) {
        snprintf(error, error_len, "%s", failure_text(failure_seen(&failed)));
        status = -1;
    }

    /* The consumer has been joined, so the writer is reachable from one thread again and
     * the run totals can be attached. */
    if (status == 0 && (cons.status < 0 ||
                        h5writer_total(p.out, OUT_UNMAPPED, unmapped) < 0)) {
        status = h5writer_fail(p.out, cfg->output_path, error, error_len);
    }

done:
    for (size_t i = 0; i < cfg->workers; i++) {
        accum_free(&workers[i].shadow);
        pairs_free(&workers[i].shadow_pairs);
        tally_scratch_destroy(workers[i].scratch);
        free(workers[i].slots);
    }

    free(workers);
    pipeline_teardown(&p);
    return status;
}
