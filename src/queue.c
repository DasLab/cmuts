/* queue.c -- bounded blocking queue, as a ring buffer under one lock.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "queue.h"

#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

struct queue {
    void          **slots;
    size_t          capacity;
    size_t          head;     /* index of the next item to leave */
    size_t          count;
    bool            closed;
    pthread_mutex_t lock;
    pthread_cond_t  not_empty;
    pthread_cond_t  not_full;
};

queue *queue_create(size_t capacity)
{
    queue *q = calloc(1, sizeof *q);
    if (!q) {
        return NULL;
    }

    q->slots = calloc(capacity, sizeof *q->slots);
    if (!q->slots) {
        free(q);
        return NULL;
    }

    q->capacity = capacity;
    pthread_mutex_init(&q->lock, NULL);
    pthread_cond_init(&q->not_empty, NULL);
    pthread_cond_init(&q->not_full, NULL);

    return q;
}

void queue_destroy(queue *q)
{
    if (!q) {
        return;
    }

    pthread_cond_destroy(&q->not_full);
    pthread_cond_destroy(&q->not_empty);
    pthread_mutex_destroy(&q->lock);
    free(q->slots);
    free(q);
}

/* ------------------------------------------------------------------------ */
/* Transfer, with the lock already held                                      */
/* ------------------------------------------------------------------------ */

static size_t fill(queue *q, void *const *items, size_t n)
{
    size_t room = q->capacity - q->count;
    size_t take = n < room ? n : room;

    for (size_t i = 0; i < take; i++) {
        q->slots[(q->head + q->count + i) % q->capacity] = items[i];
    }

    q->count += take;
    return take;
}

static size_t drain(queue *q, void **items, size_t n)
{
    size_t take = q->count < n ? q->count : n;

    for (size_t i = 0; i < take; i++) {
        items[i] = q->slots[(q->head + i) % q->capacity];
    }

    q->head   = (q->head + take) % q->capacity;
    q->count -= take;
    return take;
}

/* ------------------------------------------------------------------------ */
/* Public operations                                                         */
/* ------------------------------------------------------------------------ */

size_t queue_push(queue *q, void *const *items, size_t n)
{
    size_t pushed = 0;

    pthread_mutex_lock(&q->lock);
    while (pushed < n) {
        while (q->count == q->capacity && !q->closed) {
            pthread_cond_wait(&q->not_full, &q->lock);
        }

        if (q->closed) {
            break;
        }

        pushed += fill(q, items + pushed, n - pushed);

        /* One waiter rather than all of them: a fill adds at most a batch and a
         * consumer takes at most a batch, the same configured value on both sides, so
         * one waiter drains what one fill added.
         *
         * No wakeup is lost. A consumer holds the lock while it finds the queue empty
         * and begins to wait, and a producer holds it to fill, so either the consumer
         * is already waiting when the signal is sent or it has yet to look and will
         * find the items. */
        pthread_cond_signal(&q->not_empty);
    }
    pthread_mutex_unlock(&q->lock);

    return pushed;
}

void queue_push_all(queue *q, void *const *items, size_t n)
{
    if (queue_push(q, items, n) == n) {
        return;
    }

    fputs("cmuts: a closed queue refused work that had nowhere else to go\n", stderr);
    abort();
}

size_t queue_pop(queue *q, void **items, size_t n)
{
    pthread_mutex_lock(&q->lock);
    while (q->count == 0 && !q->closed) {
        pthread_cond_wait(&q->not_empty, &q->lock);
    }

    size_t got = drain(q, items, n);
    if (got) {
        pthread_cond_broadcast(&q->not_full);
    }
    pthread_mutex_unlock(&q->lock);

    return got;
}

size_t queue_try_pop(queue *q, void **items, size_t n)
{
    pthread_mutex_lock(&q->lock);
    size_t got = drain(q, items, n);
    if (got) {
        pthread_cond_broadcast(&q->not_full);
    }
    pthread_mutex_unlock(&q->lock);

    return got;
}

void queue_close(queue *q)
{
    pthread_mutex_lock(&q->lock);
    q->closed = true;
    pthread_cond_broadcast(&q->not_empty);
    pthread_cond_broadcast(&q->not_full);
    pthread_mutex_unlock(&q->lock);
}
