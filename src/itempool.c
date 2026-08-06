/* itempool.c -- work item recycling over a bounded free list.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "itempool.h"

#include <stdlib.h>

#include "queue.h"

struct itempool {
    workitem *storage;    /* the whole block, retained for teardown */
    size_t    capacity;
    queue    *available;  /* free list; blocking pop gives the pool its bound */
};

static void release_storage(itempool *p, size_t n)
{
    for (size_t i = 0; i < n; i++)
        bam_destroy1(p->storage[i].rec);

    free(p->storage);
}

static int stock_free_list(itempool *p)
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

itempool *itempool_create(size_t capacity)
{
    itempool *p = calloc(1, sizeof *p);
    if (!p)
        return NULL;

    p->capacity  = capacity;
    p->storage   = calloc(capacity, sizeof *p->storage);
    p->available = queue_create(capacity);
    if (!p->storage || !p->available) {
        itempool_destroy(p);
        return NULL;
    }

    for (size_t i = 0; i < capacity; i++) {
        p->storage[i].rec = bam_init1();
        if (!p->storage[i].rec) {
            release_storage(p, i);
            p->storage = NULL;
            itempool_destroy(p);
            return NULL;
        }
    }

    if (stock_free_list(p) < 0) {
        itempool_destroy(p);
        return NULL;
    }

    return p;
}

void itempool_destroy(itempool *p)
{
    if (!p)
        return;

    if (p->storage)
        release_storage(p, p->capacity);

    queue_destroy(p->available);
    free(p);
}

size_t itempool_take_many(itempool *p, void **items, size_t n)
{
    return queue_pop(p->available, items, n);
}

void itempool_give_many(itempool *p, void **items, size_t n)
{
    for (size_t i = 0; i < n; i++)
        ((workitem *)items[i])->ctx = NULL;

    queue_push_all(p->available, items, n);
}
