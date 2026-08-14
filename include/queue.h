/* queue.h -- bounded blocking queue with batched transfer.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

/* A fixed-capacity queue, safe for many producers and many consumers. The capacity is what
 * bounds the pipeline's memory: a producer that outruns its consumers blocks, and the queue
 * never grows.
 *
 * Items move in batches, so that a run of millions of records costs one lock acquisition
 * per batch and not one per item. */
typedef struct queue queue;

queue *queue_create(size_t capacity);
void   queue_destroy(queue *q);

/* Appends n items, blocking while the queue is full. Returns the number accepted, which is
 * short of n only when the queue is closed. */
size_t queue_push(queue *q, void *const *items, size_t n);

/* Appends every item, aborting the program on a short push.
 *
 * A queue rejects items only once closed, and every queue here closes strictly after its
 * producers finish, so a short push cannot happen and could not be handled if it did. A
 * dropped item is a unit of work or a pooled object, giving either an output silently
 * missing rows or a wait that never ends. */
void queue_push_all(queue *q, void *const *items, size_t n);

/* Removes up to n items, blocking until at least one arrives. Returns 0 only once the queue
 * is both closed and drained. */
size_t queue_pop(queue *q, void **items, size_t n);

/* Removes up to n items without blocking. A return of 0 means the queue was empty at the
 * moment of the call, and does not indicate whether more is coming. */
size_t queue_try_pop(queue *q, void **items, size_t n);

/* Signals the end of input, releasing every consumer once the queue drains. */
void   queue_close(queue *q);
