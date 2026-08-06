/* pipeline.c -- the loader, the worker pool, and the completion consumer.
 *
 * The file is read once, sequentially, by a single loader. Each read is copied
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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bam.h"
#include "h5writer.h"
#include "itempool.h"
#include "metadata.h"
#include "process.h"
#include "queue.h"
#include "refctx.h"
#include "refseq.h"

/* Completed references are drained in small groups. The consumer is never the
 * bottleneck, so this only saves a lock acquisition per reference. */
#define COMPLETION_BATCH 16

typedef struct {
    cm_bam_reader *bam;
    refseq_source *refs;
    queue         *work;       /* loader -> workers */
    queue         *completed;  /* last handle dropped -> consumer */
    itempool      *items;
    ctxpool       *contexts;
    h5writer      *out;
    filter_config  filter;
    size_t         batch;
    size_t         ref_cap;    /* longest reference, sizing every accumulator */
} pipeline;

static void finish_reference(pipeline *p, refctx *ctx)
{
    void *handle = ctx;

    queue_push(p->completed, &handle, 1);
}

/* ------------------------------------------------------------------------ */
/* Workers                                                                   */
/* ------------------------------------------------------------------------ */

typedef struct {
    pipeline *pipe;
    refctx   *held;       /* reference the shadow holds a handle for, or NULL */
    accum     shadow;
    void    **slots;      /* batch buffer, allocated before the thread starts */
    pthread_t thread;
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

/* Processes a run of reads that all belong to one reference, which lets the
 * whole run cost a single handle release. */
static void process_run(worker *w, void **slots, size_t n)
{
    workitem       *first = slots[0];
    refctx         *ctx   = first->ctx;
    cm_fasta_record ref;

    shadow_switch(w, ctx);
    refctx_sequence(ctx, &ref);

    for (size_t i = 0; i < n; i++) {
        const workitem *item = slots[i];
        cm_bam_record   read;

        cm_bam_record_view(item->rec, &read);
        process(&read, &ref, &w->shadow);
    }

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

static void dispatch_batch(pipeline *p, refctx *ctx, void **batch, size_t *n)
{
    if (*n == 0)
        return;

    refctx_acquire(ctx, (int)*n);
    queue_push(p->work, batch, *n);
    *n = 0;
}

/* Deposits the reads rejected for this reference before the loader's handle
 * goes, so the count cannot arrive after the reference has been written. */
static void close_reference(pipeline *p, refctx *ctx, size_t filtered)
{
    if (!ctx)
        return;

    if (filtered)
        refctx_add_scalar(ctx, ACCUM_FILTERED, (double)filtered);

    if (refctx_release(ctx, 1))
        finish_reference(p, ctx);
}

static refctx *open_reference(pipeline *p, int32_t tid)
{
    const cm_fasta_record *seq = refseq_advance(p->refs, tid);
    if (!seq)
        return NULL;

    refctx *ctx = ctxpool_take(p->contexts);
    if (!ctx)
        return NULL;

    refctx_open(ctx, tid, cm_bam_refname(p->bam, tid), seq);
    return ctx;
}

/* Leaves through a single exit, so that a partly filled batch is always
 * dispatched and the reference in hand always has the loader's handle
 * released. Otherwise a failure part way through would leave a reference that
 * can never complete. */
static int loader_main(pipeline *p, size_t *unmapped, char *error, size_t error_len)
{
    void        **batch    = calloc(p->batch, sizeof *batch);
    void        **spare    = calloc(p->batch, sizeof *spare);
    size_t        held     = 0;  /* carriers drawn from the pool, not yet used */
    refctx       *current  = NULL;
    size_t        n        = 0;
    size_t        rejected = 0;  /* filtered reads for the reference in hand */
    cm_bam_record rec;
    int           status   = CM_ITER_EOF;
    int           result   = 0;

    if (!batch || !spare) {
        snprintf(error, error_len, "out of memory");
        free(batch);
        free(spare);
        return -1;
    }

    while ((status = cm_bam_next(p->bam, &rec)) == CM_ITER_OK) {
        /* Unmapped reads align to no reference, so they are counted for the
         * run as a whole and go no further. */
        if (rec.flag & BAM_FUNMAP) {
            (*unmapped)++;
            continue;
        }

        /* The reference is opened before the filter is applied, so that reads
         * rejected here can still be counted against it: a reference whose
         * reads were all rejected is a different thing from one that received
         * none at all. */
        if (!current || rec.tid != current->tid) {
            dispatch_batch(p, current, batch, &n);
            close_reference(p, current, rejected);
            rejected = 0;

            current = open_reference(p, rec.tid);
            if (!current) {
                const char *cause = refseq_error(p->refs);

                snprintf(error, error_len, "%s", cause ? cause : "no context available");
                result = -1;
                break;
            }
        }

        if (!filter_accepts(&p->filter, &rec)) {
            rejected++;
            continue;
        }

        /* Carriers are drawn a batch at a time and spent one by one, so the
         * pool's lock is taken once per batch rather than once per read. A
         * carrier only ever moves from here into the dispatch batch, so the
         * two together never hold more than a refill's worth. */
        if (held == 0)
            held = itempool_take_many(p->items, spare, p->batch);

        if (held == 0) {
            snprintf(error, error_len, "no carrier available for a read");
            result = -1;
            break;
        }

        workitem *item = spare[--held];

        if (!bam_copy1(item->rec, cm_bam_raw(p->bam))) {
            spare[held++] = item;
            snprintf(error, error_len, "unable to copy alignment record");
            result = -1;
            break;
        }

        item->ctx  = current;
        batch[n++] = item;

        if (n == p->batch)
            dispatch_batch(p, current, batch, &n);
    }

    dispatch_batch(p, current, batch, &n);
    close_reference(p, current, rejected);

    /* Whatever the reservoir still holds goes back, on every path out. */
    if (held)
        itempool_give_many(p->items, spare, held);

    free(spare);
    free(batch);

    if (result == 0 && status == CM_ITER_ERROR) {
        snprintf(error, error_len, "%s", cm_bam_error(p->bam));
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
        .live_refs      = 0,  /* derived from the longest reference */
        .filter         = { .min_mapq = 0 },
    };
}

/* How far the loader may run ahead is what absorbs a worker that stalls on one
 * read: with too few contexts the loader blocks on the pool and the whole
 * pipeline waits behind the straggler. A context costs memory in proportion to
 * the longest reference, so the count comes from a byte budget rather than a
 * fixed number -- many short references get generous slack, while a handful of
 * very long ones cannot exhaust memory. Two is the least that still lets one
 * reference drain while the next is being loaded. */
#define TARGET_CONTEXT_BYTES (32u << 20)
#define MIN_LIVE_REFS 2
#define MAX_LIVE_REFS 256

static size_t derive_live_refs(size_t ref_cap)
{
    size_t per_context = ref_cap + 1 + accum_bytes(ref_cap);
    size_t n           = TARGET_CONTEXT_BYTES / per_context;

    if (n < MIN_LIVE_REFS)
        return MIN_LIVE_REFS;
    if (n > MAX_LIVE_REFS)
        return MAX_LIVE_REFS;

    return n;
}

static void pipeline_teardown(pipeline *p)
{
    h5writer_close(p->out);
    ctxpool_destroy(p->contexts);
    itempool_destroy(p->items);
    queue_destroy(p->completed);
    queue_destroy(p->work);
    refseq_close(p->refs);
    cm_bam_close(p->bam);
}

static int open_inputs(pipeline *p, const pipeline_config *cfg, char *error, size_t error_len)
{
    p->bam = cm_bam_open(cfg->bam_path);
    if (!p->bam) {
        snprintf(error, error_len, "%s: unable to open", cfg->bam_path);
        return -1;
    }

    if (!cm_bam_is_coordinate_sorted(p->bam)) {
        snprintf(error, error_len,
                 "%s is not coordinate sorted; references would stay live to the "
                 "end of the file", cfg->bam_path);
        return -1;
    }

    if (cm_bam_nref(p->bam) < 1) {
        snprintf(error, error_len, "%s declares no references", cfg->bam_path);
        return -1;
    }

    if (cm_bam_set_threads(p->bam, cfg->decode_threads) < 0) {
        snprintf(error, error_len, "%s", cm_bam_error(p->bam));
        return -1;
    }

    p->refs = refseq_open(cfg->fasta_path, p->bam);
    if (!p->refs) {
        snprintf(error, error_len, "%s: unable to open", cfg->fasta_path);
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
    size_t live;

    p->batch   = cfg->batch;
    p->filter  = cfg->filter;
    p->ref_cap = (size_t)cm_bam_max_reflen(p->bam);
    live       = cfg->live_refs ? cfg->live_refs : derive_live_refs(p->ref_cap);

    p->work      = queue_create(cfg->queue_capacity);
    p->completed = queue_create(live);
    p->items     = itempool_create(carriers);
    p->contexts  = ctxpool_create(live, p->ref_cap);

    if (!p->work || !p->completed || !p->items || !p->contexts) {
        snprintf(error, error_len, "out of memory building the pipeline");
        return -1;
    }

    return 0;
}

static int open_output(pipeline *p, const pipeline_config *cfg, char *error, size_t error_len)
{
    p->out = h5writer_create(cfg->output_path, cm_bam_nref(p->bam), p->ref_cap);
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
        workers[i].pipe  = p;
        workers[i].slots = calloc(p->batch, sizeof *workers[i].slots);

        if (!workers[i].slots || accum_alloc(&workers[i].shadow, p->ref_cap) < 0) {
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

    if (open_inputs(&p, cfg, error, error_len) < 0 ||
        build_buffers(&p, cfg, error, error_len) < 0 ||
        open_output(&p, cfg, error, error_len) < 0)
        goto done;

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

    /* The consumer has been joined, so the writer is again reachable from one
     * thread only and the run totals can be attached. */
    if (status == 0 && (cons.status < 0 || metadata_write_run(p.out, unmapped) < 0)) {
        snprintf(error, error_len, "%s: %s", cfg->output_path, h5writer_error(p.out));
        status = -1;
    }

done:
    for (size_t i = 0; i < cfg->workers; i++) {
        accum_free(&workers[i].shadow);
        free(workers[i].slots);
    }

    free(workers);
    pipeline_teardown(&p);
    return status;
}
