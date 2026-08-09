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
#include <sys/stat.h>
#include <unistd.h>

#include "bamstream.h"
#include "h5writer.h"
#include "itempool.h"
#include "metadata.h"
#include "progress.h"
#include "queue.h"
#include "refctx.h"
#include "refseq.h"
#include "tally.h"

/* Completed references are drained in small groups. The consumer is never the
 * bottleneck, so this only saves a lock acquisition per reference. */
#define COMPLETION_BATCH 16

typedef struct {
    cm_bam_stream *bam;
    refseq_source *refs;
    queue         *work;       /* loader -> workers */
    queue         *completed;  /* last handle dropped -> consumer */
    itempool      *items;
    ctxpool       *contexts;
    h5writer      *out;
    progress      *bar;
    tally_tables   tally_tables;
    /* The first way a worker found to stop, which every worker writes and the
     * loader reads: one failure is enough, and which of them arrived first is
     * all the run has to report. */
    _Atomic int    failure;   /* a phmm_status */
    bool           may_replace;
    filter_config  filter_config;
    size_t         batch;
    size_t         ref_cap;    /* longest reference, sizing every accumulator */
} pipeline;

static void finish_reference(pipeline *p, refctx *ctx)
{
    void *handle = ctx;

    queue_push_all(p->completed, &handle, 1);
}

/* ------------------------------------------------------------------------ */
/* Workers                                                                   */
/* ------------------------------------------------------------------------ */

typedef struct {
    pipeline      *pipe;
    refctx        *held;    /* reference the shadow holds a handle for, or NULL */
    accum          shadow;
    tally_scratch *scratch; /* what the processing step works in */
    void         **slots;   /* batch buffer, allocated before the thread starts */
    pthread_t      thread;
} worker;

/* Merges the shadow into its reference and releases the handle that protected
 * it. The shadow is left entirely zero: it is cleared over exactly the extent
 * it was dirtied, and was zero elsewhere by the same argument on the previous
 * flush. */
static void shadow_flush(worker *w)
{
    refctx *ctx = w->held;

    if (!ctx)
        return;

    refctx_merge(ctx, &w->shadow);
    accum_zero(&w->shadow, ctx->len);
    w->held = NULL;

    if (refctx_release(ctx, 1))
        finish_reference(w->pipe, ctx);
}

static void shadow_switch(worker *w, refctx *ctx)
{
    if (w->held == ctx)
        return;

    shadow_flush(w);
    refctx_acquire(ctx, 1);
    w->held = ctx;
}

/* Every read of a run, into the worker's own shadow. A read the tally cannot
 * count is one no later read would fare better on, so what it found is put
 * where the loader will see it and the batch is seen through regardless. */
static void count_run(worker *w, void **slots, size_t n,
                      const cm_fasta_record *ref)
{
    for (size_t i = 0; i < n; i++) {
        const workitem *item = slots[i];
        cm_bam_record   read;
        phmm_status     status;

        cm_bam_record_view(item->rec, &read);
        status = tally(&read, ref, &w->pipe->tally_tables, w->scratch,
                       &w->shadow);

        if (status != PHMM_OK)
            atomic_store(&w->pipe->failure, status);
    }
}

/* Processes a run of reads that all belong to one reference, which lets the
 * whole run cost a single handle release. */
static void process_run(worker *w, void **slots, size_t n)
{
    workitem       *first = slots[0];
    refctx         *ctx   = first->ctx;
    cm_fasta_record ref;

    shadow_switch(w, ctx);
    refctx_sequence(ctx, &ref);

    /* A run the failure reached first is unwound and not counted: nothing will
     * read what it would have contributed, and the handles it holds are owed
     * back whether it contributes or not. */
    if (atomic_load(&w->pipe->failure) == PHMM_OK)
        count_run(w, slots, n, &ref);

    /* The whole run is finished at the same moment, so it goes back in one
     * piece. Nothing reads these carriers afterwards: process_batch scans for
     * runs only at or beyond the current head, and every item before it has
     * already been processed. */
    itempool_give_many(w->pipe->items, slots, n);

    /* Cannot be the last handle: the shadow holds one on this reference. */
    if (refctx_release(ctx, (int)n))
        finish_reference(w->pipe, ctx);
}

static void process_batch(worker *w, void **slots, size_t n)
{
    for (size_t i = 0; i < n; ) {
        const workitem *head = slots[i];
        size_t          run  = 1;

        while (i + run < n && ((const workitem *)slots[i + run])->ctx == head->ctx)
            run++;

        process_run(w, slots + i, run);
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
            /* Flush before blocking: an idle worker holding a shadow would pin
             * its reference indefinitely, and the loader would eventually
             * stall waiting for a context to come free. */
            shadow_flush(w);

            n = queue_pop(w->pipe->work, w->slots, batch);
            if (n == 0)
                break;
        }

        process_batch(w, w->slots, n);
    }

    shadow_flush(w);
    return NULL;
}

/* ------------------------------------------------------------------------ */
/* Loader                                                                    */
/* ------------------------------------------------------------------------ */

/* What the loader has in hand.
 *
 * Each of the three is a pair of values that only make sense together: a
 * reference and the reads it has turned away, a batch bound for the workers
 * and how much of it is filled, a reservoir of carriers and how many remain.
 * Kept as locals they were six variables to hold in step by hand. */
typedef struct {
    pipeline *pipe;

    refctx *reference;  /* the one being filled, or none yet */
    size_t  rejected;   /* its reads the filter turned away */

    void  **batch;      /* reads bound for the workers */
    size_t  queued;

    void  **spare;      /* carriers drawn from the pool, not yet spent */
    size_t  held;
} loader;

static int loader_open(loader *l, pipeline *p)
{
    *l = (loader){ .pipe = p };

    l->batch = calloc(p->batch, sizeof *l->batch);
    l->spare = calloc(p->batch, sizeof *l->spare);

    if (l->batch && l->spare)
        return 0;

    free(l->batch);
    free(l->spare);
    return -1;
}

/* Hands what is queued to the workers, taking a handle on the reference for
 * each read so that it cannot be finished while any of them is in transit. */
static void loader_dispatch(loader *l)
{
    if (l->queued == 0)
        return;

    refctx_acquire(l->reference, (int)l->queued);
    queue_push_all(l->pipe->work, l->batch, l->queued);
    l->queued = 0;
}

/* Lets go of the reference in hand: what is queued for it goes first, then the
 * count of what it turned away, and the loader's own handle last. The count
 * has to arrive before the handle goes, or the reference could be written
 * without it. */
static void loader_leave_reference(loader *l)
{
    refctx *ctx = l->reference;

    if (!ctx)
        return;

    loader_dispatch(l);

    if (l->rejected)
        refctx_add_scalar(ctx, ACCUM_FILTERED, (double)l->rejected);

    l->reference = NULL;
    l->rejected  = 0;

    if (refctx_release(ctx, 1))
        finish_reference(l->pipe, ctx);
}

static refctx *open_reference(pipeline *p, int32_t tid)
{
    const cm_fasta_record *seq = refseq_advance(p->refs, tid);
    if (!seq)
        return NULL;

    refctx *ctx = ctxpool_take(p->contexts);
    if (!ctx)
        return NULL;

    refctx_open(ctx, tid, cm_bam_stream_refname(p->bam, tid), seq);
    return ctx;
}

/* Moves to the reference a read belongs to, unless it is already the one in
 * hand. This happens before the filter is applied, so that a read rejected
 * there can still be counted against its reference: one whose reads were all
 * turned away is a different thing from one that received none at all. */
static bool loader_on_reference(loader *l, int32_t tid)
{
    if (l->reference && l->reference->tid == tid)
        return true;

    loader_leave_reference(l);
    l->reference = open_reference(l->pipe, tid);

    return l->reference != NULL;
}

/* Carriers are drawn a batch at a time and spent one by one, so the pool's
 * lock is taken once per batch rather than once per read. One only ever moves
 * from the reservoir into the batch, so the two together never hold more than
 * a refill's worth. */
static workitem *loader_carrier(loader *l)
{
    if (l->held == 0)
        l->held = itempool_take_many(l->pipe->items, l->spare, l->pipe->batch);

    return l->held ? l->spare[--l->held] : NULL;
}

/* Takes a copy of the record just read, since the reader overwrites its own on
 * the next advance, and queues it for the workers. */
static int loader_take_read(loader *l)
{
    workitem *item = loader_carrier(l);

    if (!item)
        return -1;

    if (!bam_copy1(item->rec, cm_bam_stream_raw(l->pipe->bam))) {
        l->spare[l->held++] = item;
        return -1;
    }

    item->ctx             = l->reference;
    l->batch[l->queued++] = item;

    if (l->queued == l->pipe->batch)
        loader_dispatch(l);

    return 0;
}

/* Everything in hand goes back, whichever way the loop was left. */
static void loader_finish(loader *l)
{
    loader_leave_reference(l);

    if (l->held)
        itempool_give_many(l->pipe->items, l->spare, l->held);

    free(l->spare);
    free(l->batch);
}

/* Reads the file once, in order, dispatching what survives. */
static int loader_main(pipeline *p, size_t *unmapped, char *error, size_t error_len)
{
    loader        l;
    cm_bam_record rec;
    int           status = CM_ITER_EOF;
    int           result = 0;

    if (loader_open(&l, p) < 0) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    while ((status = cm_bam_stream_next(p->bam, &rec)) == CM_ITER_OK) {
        progress_follow(p->bar);

        /* A worker has found something no later read will mend, so there is
         * nothing to be had from reading them. Stopping is not failing: what
         * went wrong is the worker's to report, not the loader's. */
        if (atomic_load(&p->failure) != PHMM_OK)
            break;

        /* Unmapped reads align to no reference, so they are counted for the
         * run as a whole and go no further. */
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
    pipeline *pipe;
    int       status;  /* first write failure, if any */
    pthread_t thread;
} consumer;

static void *consumer_main(void *arg)
{
    consumer *c = arg;
    void     *slots[COMPLETION_BATCH];
    size_t    n;

    while ((n = queue_pop(c->pipe->completed, slots, COMPLETION_BATCH)) > 0)
        for (size_t i = 0; i < n; i++) {
            refctx *ctx = slots[i];

            /* Keep draining after a failure: the loader and workers must not
             * be left blocked on a queue nobody is emptying. */
            if (c->status == 0 &&
                h5writer_row(c->pipe->out, ctx->tid, ctx->len, &ctx->acc) < 0)
                c->status = -1;

            ctxpool_give(c->pipe->contexts, ctx);
        }

    return NULL;
}

/* ------------------------------------------------------------------------ */
/* Assembly                                                                  */
/* ------------------------------------------------------------------------ */

pipeline_config pipeline_defaults(void)
{
    return (pipeline_config){
        .workers        = 4,
        .decode_threads = 4,
        .queue_capacity = 4096,
        .batch          = 64,
        .live_refs      = 64,
        .filter_config  = filter_defaults(),
        .tally_config   = tally_defaults(),
    };
}

static void pipeline_teardown(pipeline *p)
{
    progress_finish(p->bar);
    h5writer_close(p->out);
    ctxpool_destroy(p->contexts);
    itempool_destroy(p->items);
    queue_destroy(p->completed);
    queue_destroy(p->work);
    refseq_close(p->refs);
    cm_bam_stream_close(p->bam);
}

/* What is worth refusing is a file with something in it. A path that exists
 * but is empty is what mktemp and shell redirection leave behind, and there is
 * nothing there to lose. */
static bool holds_data(const char *path)
{
    struct stat info;

    return stat(path, &info) == 0 && info.st_size > 0;
}

/* Checked before anything is opened, so a mistyped path costs nothing and no
 * previous result is at risk while the inputs are still being read. */
static int check_output(pipeline *p, const pipeline_config *cfg,
                        char *error, size_t error_len)
{
    if (!cfg->overwrite && holds_data(cfg->output_path)) {
        snprintf(error, error_len,
                 "%s already holds data; pass --overwrite to replace it",
                 cfg->output_path);
        return -1;
    }

    /* Where nothing is at the path, the create stays exclusive, so a file
     * appearing in between is not quietly replaced. Where an empty one is
     * already there it has to be truncated instead. */
    p->may_replace = cfg->overwrite || access(cfg->output_path, F_OK) == 0;
    return 0;
}

static int open_inputs(pipeline *p, const pipeline_config *cfg, char *error, size_t error_len)
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

    p->refs = refseq_open(cfg->fasta_path, p->bam, &why);
    if (!p->refs) {
        snprintf(error, error_len, "%s: %s", cfg->fasta_path, why);
        return -1;
    }

    return 0;
}

static int build_buffers(pipeline *p, const pipeline_config *cfg, char *error, size_t error_len)
{
    /* Enough carriers for a full queue, a batch in the hands of every worker,
     * and two for the loader: one being filled and one held in reserve, since
     * a short refill can leave it holding part of each. */
    size_t carriers = cfg->queue_capacity + (cfg->workers + 2) * cfg->batch;

    p->batch         = cfg->batch;
    p->filter_config = cfg->filter_config;
    p->ref_cap       = (size_t)cm_bam_stream_max_reflen(p->bam);

    p->work      = queue_create(cfg->queue_capacity);
    p->completed = queue_create(cfg->live_refs);
    p->items     = itempool_create(carriers);
    p->contexts  = ctxpool_create(cfg->live_refs, p->ref_cap);

    if (!p->work || !p->completed || !p->items || !p->contexts) {
        snprintf(error, error_len, "out of memory building the pipeline");
        return -1;
    }

    return 0;
}

static int open_output(pipeline *p, const pipeline_config *cfg, char *error, size_t error_len)
{
    p->out = h5writer_create(cfg->output_path, cm_bam_stream_nref(p->bam), p->ref_cap,
                             p->may_replace);
    if (!p->out) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    if (h5writer_error(p->out) || metadata_write_names(p->out, p->bam) < 0) {
        const char *cause = h5writer_error(p->out);

        snprintf(error, error_len, "%s: %s", cfg->output_path,
                 cause ? cause : "unable to write the reference names");
        return -1;
    }

    return 0;
}

/* Starts as many workers as it can, reporting through started how many are
 * running so that the caller joins exactly those. Every allocation a worker
 * needs is made here rather than inside the thread, so that running out of
 * memory is reported instead of leaving a thread that exits at once and a
 * loader that waits forever for it. */
static int start_workers(worker *workers, size_t n, pipeline *p, size_t *started,
                         char *error, size_t error_len)
{
    *started = 0;

    for (size_t i = 0; i < n; i++) {
        workers[i].pipe    = p;
        workers[i].slots   = calloc(p->batch, sizeof *workers[i].slots);
        workers[i].scratch = tally_scratch_create();

        if (!workers[i].slots || !workers[i].scratch
            || accum_alloc(&workers[i].shadow, p->ref_cap) < 0) {
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

int pipeline_run(const pipeline_config *cfg, char *error, size_t error_len)
{
    pipeline  p        = { 0 };
    consumer  cons     = { 0 };
    size_t    unmapped = 0;
    worker   *workers  = calloc(cfg->workers, sizeof *workers);
    size_t    started  = 0;
    int       status   = -1;

    if (!workers) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    if (check_output(&p, cfg, error, error_len) < 0 ||
        open_inputs(&p, cfg, error, error_len) < 0 ||
        build_buffers(&p, cfg, error, error_len) < 0 ||
        open_output(&p, cfg, error, error_len) < 0)
        goto done;

    tally_tables_build(&p.tally_tables, &cfg->tally_config);

    /* Once nothing is left that could fail before a read is taken. */
    p.bar = progress_start(p.bam);

    cons.pipe = &p;
    if (pthread_create(&cons.thread, NULL, consumer_main, &cons) != 0) {
        snprintf(error, error_len, "unable to start the consumer");
        goto done;
    }

    if (start_workers(workers, cfg->workers, &p, &started, error, error_len) == 0)
        status = loader_main(&p, &unmapped, error, error_len);

    queue_close(p.work);
    for (size_t i = 0; i < started; i++)
        pthread_join(workers[i].thread, NULL);

    queue_close(p.completed);
    pthread_join(cons.thread, NULL);

    /* Asked once every worker has stopped, so that a failure on the last batch
     * is caught as surely as one the loader saw in time to stop for. A loader
     * that failed on its own account has already said something more particular
     * than this could. */
    if (status == 0 && atomic_load(&p.failure) != PHMM_OK) {
        snprintf(error, error_len, "%s",
                 atomic_load(&p.failure) == PHMM_NO_MEMORY
                     ? "out of memory marginalizing a read"
                     : "a marginalization did not hold together");
        status = -1;
    }

    /* The consumer has been joined, so the writer is again reachable from one
     * thread only and the run totals can be attached. */
    if (status == 0 && (cons.status < 0 || metadata_write_run(p.out, unmapped) < 0)) {
        snprintf(error, error_len, "%s: %s", cfg->output_path, h5writer_error(p.out));
        status = -1;
    }

done:
    for (size_t i = 0; i < cfg->workers; i++) {
        accum_free(&workers[i].shadow);
        tally_scratch_destroy(workers[i].scratch);
        free(workers[i].slots);
    }

    free(workers);
    pipeline_teardown(&p);
    return status;
}
