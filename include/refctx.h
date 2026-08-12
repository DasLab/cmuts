/* refctx.h -- per-reference accumulation state and its completion protocol.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <pthread.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "accum.h"
#include "fasta.h"

/* Everything a reference accumulates while its reads are in flight.
 *
 * Completion is tracked with a handle count rather than a done flag: the loader holds one
 * while it is still dispatching, each read in transit holds one, and each worker holding
 * unmerged contributions in a shadow accumulator holds one. The reference is finished
 * exactly when the count falls to zero, which leaves no ordering question between the
 * loader closing a reference and the workers draining it. */
typedef struct refctx {
    int32_t         tid;   /* also the output row index */
    const char     *name;  /* borrowed from the BAM header, stable for the run */
    char           *seq;   /* owned; the FASTA reader's buffer does not persist */
    size_t          len;
    accum           acc;
    atomic_int      handles;
    pthread_mutex_t lock;  /* guards acc */
} refctx;

void refctx_acquire(refctx *ctx, int n);

/* Drops n handles, returning true to whichever caller drops the last one. That caller then
 * owns the reference and must pass it on to be written. */
bool refctx_release(refctx *ctx, int n);

/* ctx->acc += src, serialized against the other workers. The lock is held for the merge
 * and not for the processing, which happens outside it. */
void refctx_merge(refctx *ctx, const accum *src);

/* Adds to one scalar field, for counts arising outside the processing step and so having
 * no accumulator of their own to merge from. */
void refctx_add_scalar(refctx *ctx, accum_field_id id, double value);

/* Presents the reference sequence in the form the processing step expects. */
void refctx_sequence(const refctx *ctx, cm_fasta_record *out);

/* A fixed set of contexts. Its size caps how many references may be in flight, bounding
 * memory independently of how many references the file declares.
 *
 * A pooled context always has a zeroed accumulator, an invariant ctxpool_give restores, so
 * a context taken from the pool is ready to accumulate whatever the previous reference
 * left behind. */
typedef struct ctxpool ctxpool;

/* ref_cap is the longest reference in the file. Every context is sized to it so that any
 * two accumulators share a layout and may be merged. */
ctxpool *ctxpool_create(size_t capacity, size_t ref_cap);
void     ctxpool_destroy(ctxpool *p);

/* Takes a context, blocking while every one is live. */
refctx *ctxpool_take(ctxpool *p);
void    ctxpool_give(ctxpool *p, refctx *ctx);

/* Readies a context taken from the pool, with the loader's handle held. */
void refctx_open(refctx *ctx, int32_t tid, const char *name, const cm_fasta_record *seq);
