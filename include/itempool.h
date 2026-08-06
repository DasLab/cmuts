/* itempool.h -- recycled carriers for reads in transit.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include <htslib/sam.h>

struct refctx;

/* One read on its way to a worker: the reference it belongs to, and a private
 * copy of the alignment, since the reader overwrites its own record on every
 * advance. The copy is made with bam_copy1 into a record that the pool keeps
 * alive for the whole run, so a file of millions of reads performs no
 * per-read allocation. */
typedef struct {
    struct refctx *ctx;
    bam1_t        *rec;
} workitem;

typedef struct itempool itempool;

itempool *itempool_create(size_t capacity);
void      itempool_destroy(itempool *p);

/* Carriers move in groups for the same reason reads do: taken one at a time,
 * the free list's lock is acquired twice per read, which costs more than the
 * processing it protects and caps the pipeline well short of its throughput.
 *
 * Takes up to n, blocking only while the pool is completely exhausted, and
 * returns how many were obtained. That may be fewer than asked for, which is
 * not an error; it is what bounds the reads resident in memory. */
size_t itempool_take_many(itempool *p, void **items, size_t n);
void   itempool_give_many(itempool *p, void **items, size_t n);
