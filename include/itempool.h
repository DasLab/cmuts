/* itempool.h -- recycled carriers for reads in transit.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include <htslib/sam.h>

struct refctx;
struct h5chunk;

/* One unit of work on its way to a worker. Almost always a read: the reference it belongs
 * to, and a private copy of the alignment, the reader overwriting its own record on every
 * advance. The copy is made with bam_copy1 into a record the pool keeps alive for the
 * whole run, so a file of millions of reads performs no per-read allocation.
 *
 * An item with chunk set instead carries an output chunk to filter for the writing
 * thread. Such an item never comes from the pool: it is allocated for its one trip and
 * freed by the worker that takes it. */
typedef struct {
    struct refctx  *ctx;
    bam1_t         *rec;
    struct h5chunk *chunk;
} workitem;

typedef struct itempool itempool;

itempool *itempool_create(size_t capacity);
void      itempool_destroy(itempool *p);

/* Takes up to n carriers, blocking only while the pool is entirely exhausted, and returns
 * how many were obtained. Fewer than asked for is not an error; it is what bounds the reads
 * resident in memory.
 *
 * Taken in groups for the same reason reads are, so that the free list's lock is acquired
 * once per group and not twice per read. */
size_t itempool_take_many(itempool *p, void **items, size_t n);
void   itempool_give_many(itempool *p, void **items, size_t n);
