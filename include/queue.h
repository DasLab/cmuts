/* queue.h -- bounded blocking queue with batched transfer.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

/* A fixed-capacity queue, safe for many producers and many consumers. The
 * capacity is what bounds the pipeline's memory: a producer that outruns its
 * consumers blocks rather than growing the queue.
 *
 * Items move in batches so that a run of millions of records costs one lock
 * acquisition per batch rather than one per item. */
typedef struct queue queue;

queue *queue_create(size_t capacity);
void   queue_destroy(queue *q);

/* Appends n items, blocking while the queue is full. Returns the number
 * accepted, which is short of n only when the queue is closed. */
size_t queue_push(queue *q, void *const *items, size_t n);

/* Removes up to n items, blocking until at least one arrives. Returns 0 only
 * once the queue is both closed and drained. */
size_t queue_pop(queue *q, void **items, size_t n);

/* Removes up to n items without ever blocking. A return of 0 means the queue
 * was empty at the moment of the call, and says nothing about whether more is
 * coming. */
size_t queue_try_pop(queue *q, void **items, size_t n);

/* Signals the end of input, releasing every consumer once the queue drains. */
void   queue_close(queue *q);
