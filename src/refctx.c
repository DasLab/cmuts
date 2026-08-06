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
    queue  *available;  /* free list; blocking pop gives the pool its bound */
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

void refctx_merge(refctx *ctx, const accum *src)
{
    pthread_mutex_lock(&ctx->lock);
    accum_add(&ctx->acc, src, ctx->len);
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

static void release_storage(ctxpool *p, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        pthread_mutex_destroy(&p->storage[i].lock);
        accum_free(&p->storage[i].acc);
        free(p->storage[i].seq);
    }

    free(p->storage);
}

static int build_context(refctx *ctx, size_t ref_cap)
{
    ctx->seq = malloc(ref_cap + 1);
    if (!ctx->seq)
        return -1;

    if (accum_alloc(&ctx->acc, ref_cap) < 0) {
        free(ctx->seq);
        ctx->seq = NULL;
        return -1;
    }

    pthread_mutex_init(&ctx->lock, NULL);
    atomic_init(&ctx->handles, 0);
    return 0;
}

static int stock_free_list(ctxpool *p)
{
    void **handles = calloc(p->capacity, sizeof *handles);
    if (!handles)
        return -1;

    for (size_t i = 0; i < p->capacity; i++)
        handles[i] = &p->storage[i];

    size_t pushed = queue_push(p->available, handles, p->capacity);
    free(handles);

    return pushed == p->capacity ? 0 : -1;
}

ctxpool *ctxpool_create(size_t capacity, size_t ref_cap)
{
    ctxpool *p = calloc(1, sizeof *p);
    if (!p)
        return NULL;

    p->capacity  = capacity;
    p->ref_cap   = ref_cap;
    p->storage   = calloc(capacity, sizeof *p->storage);
    p->available = queue_create(capacity);
    if (!p->storage || !p->available) {
        ctxpool_destroy(p);
        return NULL;
    }

    for (size_t i = 0; i < capacity; i++) {
        if (build_context(&p->storage[i], ref_cap) < 0) {
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
    if (!p)
        return;

    if (p->storage)
        release_storage(p, p->capacity);

    queue_destroy(p->available);
    free(p);
}

refctx *ctxpool_take(ctxpool *p)
{
    refctx *ctx = NULL;

    return queue_pop(p->available, (void **)&ctx, 1) == 1 ? ctx : NULL;
}

void ctxpool_give(ctxpool *p, refctx *ctx)
{
    accum_zero(&ctx->acc, ctx->len);
    ctx->len = 0;
    queue_push(p->available, (void *const *)&ctx, 1);
}
