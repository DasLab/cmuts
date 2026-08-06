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

/* Blocks while every item is checked out. This is what bounds the number of
 * reads resident in memory, whatever the size of the file. */
workitem *itempool_take(itempool *p);
void      itempool_give(itempool *p, workitem *item);
