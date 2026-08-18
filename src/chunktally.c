/* chunktally.c -- counting references per output chunk.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "chunktally.h"

#include <pthread.h>
#include <stdbool.h>
#include <stdlib.h>

struct chunktally {
    pthread_mutex_t lock;
    size_t          rows;         /* references per chunk */
    size_t          count;        /* chunks in the output */
    int32_t        *outstanding;  /* opened here, not yet written */
    bool           *touched;      /* received at least one reference */
    bool           *handed_back;
    int64_t         opened_high;  /* the furthest chunk the loader has reached */
    bool            spent;        /* the loader will open nothing further */
    int64_t        *settled;      /* queue of chunks ready to be written */
    size_t          head;
    size_t          tail;
};

static size_t chunk_count(int32_t n_refs, size_t rows)
{
    return ((size_t)n_refs + rows - 1) / rows;
}

chunktally *chunktally_create(int32_t n_refs, size_t rows_per_chunk)
{
    chunktally *t = calloc(1, sizeof *t);

    if (!t) {
        return NULL;
    }

    t->rows        = rows_per_chunk;
    t->count       = chunk_count(n_refs, rows_per_chunk);
    t->opened_high = -1;
    t->outstanding = calloc(t->count, sizeof *t->outstanding);
    t->touched     = calloc(t->count, sizeof *t->touched);
    t->handed_back = calloc(t->count, sizeof *t->handed_back);
    t->settled     = calloc(t->count, sizeof *t->settled);

    if (!t->outstanding || !t->touched || !t->handed_back || !t->settled) {
        chunktally_destroy(t);
        return NULL;
    }

    pthread_mutex_init(&t->lock, NULL);
    return t;
}

void chunktally_destroy(chunktally *t)
{
    if (!t) {
        return;
    }

    pthread_mutex_destroy(&t->lock);
    free(t->settled);
    free(t->handed_back);
    free(t->touched);
    free(t->outstanding);
    free(t);
}

/* ------------------------------------------------------------------------ */
/* Settling, with the lock already held                                      */
/* ------------------------------------------------------------------------ */

static int64_t chunk_of(const chunktally *t, int32_t tid)
{
    return (int64_t)tid / (int64_t)t->rows;
}

/* Sealed once the loader has opened a reference beyond it, or given up on the
 * file entirely. */
static bool sealed(const chunktally *t, int64_t chunk)
{
    return t->spent || t->opened_high > chunk;
}

static void hand_back(chunktally *t, int64_t chunk)
{
    if (t->handed_back[chunk] || !t->touched[chunk] ||
        t->outstanding[chunk] != 0 || !sealed(t, chunk)) {
        return;
    }

    t->handed_back[chunk] = true;
    t->settled[t->tail++] = chunk;
}

/* Sealing a chunk can settle it without anything being written, so the range
 * the loader has just passed is offered as well as the chunk itself. */
static void hand_back_range(chunktally *t, int64_t from, int64_t through)
{
    for (int64_t c = from; c <= through && c < (int64_t)t->count; c++) {
        if (c >= 0) {
            hand_back(t, c);
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Reporting                                                                 */
/* ------------------------------------------------------------------------ */

void chunktally_expect(chunktally *t, int32_t tid)
{
    int64_t chunk = chunk_of(t, tid);

    pthread_mutex_lock(&t->lock);

    t->outstanding[chunk]++;
    t->touched[chunk] = true;

    if (chunk > t->opened_high) {
        int64_t previous = t->opened_high;

        t->opened_high = chunk;
        hand_back_range(t, previous, chunk - 1);
    }

    pthread_mutex_unlock(&t->lock);
}

void chunktally_wrote(chunktally *t, int32_t tid)
{
    int64_t chunk = chunk_of(t, tid);

    pthread_mutex_lock(&t->lock);

    if (--t->outstanding[chunk] == 0) {
        hand_back(t, chunk);
    }

    pthread_mutex_unlock(&t->lock);
}

void chunktally_no_more(chunktally *t)
{
    pthread_mutex_lock(&t->lock);

    t->spent = true;
    hand_back_range(t, 0, (int64_t)t->count - 1);

    pthread_mutex_unlock(&t->lock);
}

int64_t chunktally_take_settled(chunktally *t)
{
    int64_t chunk = -1;

    pthread_mutex_lock(&t->lock);

    if (t->head < t->tail) {
        chunk = t->settled[t->head++];
    }

    pthread_mutex_unlock(&t->lock);

    return chunk;
}
