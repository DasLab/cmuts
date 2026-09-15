/* refctx.c -- reference lifetime, merging and recycling.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "refctx.h"

#include <stdlib.h>
#include <string.h>

#include "queue.h"

struct ctxpool {
    refctx *storage;    /* the whole block, retained for teardown */
    size_t  capacity;
    size_t  ref_cap;
    queue  *available;  /* free list; a blocking pop is what bounds the pool */
};

/* ------------------------------------------------------------------------ */
/* Handles                                                                   */
/* ------------------------------------------------------------------------ */

void refctx_acquire(refctx *ctx, int n)
{
    atomic_fetch_add_explicit(&ctx->handles, n, memory_order_acq_rel);
}

bool refctx_release(refctx *ctx, int n)
{
    return atomic_fetch_sub_explicit(&ctx->handles, n, memory_order_acq_rel) == n;
}

/* ------------------------------------------------------------------------ */
/* Contents                                                                  */
/* ------------------------------------------------------------------------ */

bool refctx_accumulated(const refctx *ctx)
{
    return ctx->accumulated;
}

void refctx_merge(refctx *ctx, const accum *src, const pairs *src_pairs)
{
    pthread_mutex_lock(&ctx->lock);
    ctx->accumulated = true;
    accum_add(&ctx->acc, src, ctx->len);

    if (src_pairs) {
        pairs_add(&ctx->pr, src_pairs, ctx->len);
    }

    pthread_mutex_unlock(&ctx->lock);
}

void refctx_add_scalar(refctx *ctx, accum_field_id id, double value)
{
    pthread_mutex_lock(&ctx->lock);
    ctx->accumulated = true;
    *accum_data(&ctx->acc, id) += value;
    pthread_mutex_unlock(&ctx->lock);
}

void refctx_sequence(const refctx *ctx, cm_fasta_record *out)
{
    out->name    = ctx->name;
    out->comment = NULL;
    out->seq     = ctx->seq;
    out->len     = ctx->len;
}

void refctx_open(refctx *ctx, int32_t tid, const char *name, const cm_fasta_record *seq)
{
    ctx->tid  = tid;
    ctx->name = name;
    ctx->len  = seq->len;

    memcpy(ctx->seq, seq->seq, seq->len);
    ctx->seq[seq->len] = '\0';

    atomic_store_explicit(&ctx->handles, 1, memory_order_release);
}

/* ------------------------------------------------------------------------ */
/* Pool                                                                      */
/* ------------------------------------------------------------------------ */

/* Frees a context's buffers. Each free accepts a buffer never allocated, so this also
 * serves a context built only in part. */
static void free_buffers(refctx *ctx)
{
    pairs_free(&ctx->pr);
    accum_free(&ctx->acc);
    free(ctx->rate_storage);
    free(ctx->seq);
    ctx->rate_storage = NULL;
    ctx->seq          = NULL;
}

static void release_storage(ctxpool *p, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        pthread_mutex_destroy(&p->storage[i].lock);
        free_buffers(&p->storage[i]);
    }

    free(p->storage);
}

static int build_context(refctx *ctx, size_t ref_cap, bool pairwise)
{
    ctx->seq          = malloc(ref_cap + 1);
    ctx->rate_storage = malloc(PHMM_RATE_ARRAYS * ref_cap * sizeof *ctx->rate_storage);

    if (!ctx->seq || !ctx->rate_storage || accum_alloc(&ctx->acc, ref_cap) < 0
        || (pairwise && pairs_alloc(&ctx->pr, ref_cap) < 0)) {
        free_buffers(ctx);
        return -1;
    }

    pthread_mutex_init(&ctx->lock, NULL);
    atomic_init(&ctx->handles, 0);
    return 0;
}

static int stock_free_list(ctxpool *p)
{
    void **handles = calloc(p->capacity, sizeof *handles);
    if (!handles) {
        return -1;
    }

    for (size_t i = 0; i < p->capacity; i++) {
        handles[i] = &p->storage[i];
    }

    size_t pushed = queue_push(p->available, handles, p->capacity);
    free(handles);

    return pushed == p->capacity ? 0 : -1;
}

ctxpool *ctxpool_create(size_t capacity, size_t ref_cap, bool pairwise)
{
    ctxpool *p = calloc(1, sizeof *p);
    if (!p) {
        return NULL;
    }

    p->capacity  = capacity;
    p->ref_cap   = ref_cap;
    p->storage   = calloc(capacity, sizeof *p->storage);
    p->available = queue_create(capacity);
    if (!p->storage || !p->available) {
        ctxpool_destroy(p);
        return NULL;
    }

    for (size_t i = 0; i < capacity; i++) {
        if (build_context(&p->storage[i], ref_cap, pairwise) < 0) {
            release_storage(p, i);
            p->storage = NULL;
            ctxpool_destroy(p);
            return NULL;
        }
    }

    if (stock_free_list(p) < 0) {
        ctxpool_destroy(p);
        return NULL;
    }

    return p;
}

void ctxpool_destroy(ctxpool *p)
{
    if (!p) {
        return;
    }

    if (p->storage) {
        release_storage(p, p->capacity);
    }

    queue_destroy(p->available);
    free(p);
}

/* Takes a context from the pool, blocking while every one is live.
 *
 * The queue carries a void *, converted through a slot of its own and not by reading
 * a refctx * as though it were one. Only void * and char * are guaranteed the
 * representation of every other object pointer, so the two may not be aliased. */
refctx *ctxpool_take(ctxpool *p)
{
    void *slot = NULL;

    return queue_pop(p->available, &slot, 1) == 1 ? slot : NULL;
}

/* A pooled context always has a zeroed accumulator, so one that never accumulated is
 * already as the pool requires and is returned untouched. */
void ctxpool_give(ctxpool *p, refctx *ctx)
{
    void *slot = ctx;

    if (ctx->accumulated) {
        accum_zero(&ctx->acc, ctx->len);

        if (ctx->pr.cells) {
            pairs_zero(&ctx->pr, ctx->len);
        }

        ctx->accumulated = false;
    }

    ctx->len = 0;
    queue_push_all(p->available, &slot, 1);
}
