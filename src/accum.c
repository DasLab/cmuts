/* accum.c -- accumulator storage, laid out from the field table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "accum.h"

#include <stdlib.h>
#include <string.h>

/* Coverage and span differ over exactly what a read did not read: a deleted
 * position is spanned but not covered, and a poorly read one is spanned whole
 * but covered only in proportion to the confidence in its base. A rate over
 * coverage gives the fraction of bases actually read that disagreed; a rate
 * over span gives the fraction of reads reaching a position that found anything
 * there, which is the denominator a deletion belongs over. */
const accum_field ACCUM_FIELDS[ACCUM_N_FIELDS] = {
    [ACCUM_COVERAGE]  = { "coverage",       ACCUM_PER_BASE },
    [ACCUM_SPANNED]   = { "spanned",        ACCUM_PER_BASE },
    [ACCUM_MUTATIONS] = { "mutations",      ACCUM_PER_BASE },
    [ACCUM_LENGTHS]   = { "read_lengths",   ACCUM_PER_LENGTH },
    [ACCUM_READS]     = { "reads",          ACCUM_SCALAR   },
    [ACCUM_FILTERED]  = { "reads_filtered", ACCUM_SCALAR   },
};

/* With len == cap this also gives the field's stride in the arena. */
size_t accum_extent(accum_field_id id, size_t len)
{
    switch (ACCUM_FIELDS[id].kind) {
        case ACCUM_PER_BASE:   return len;
        case ACCUM_PER_LENGTH: return ACCUM_LENGTH_BINS(len);
        default:               return 1;
    }
}

static size_t arena_extent(size_t len)
{
    size_t total = 0;

    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++)
        total += accum_extent(id, len);

    return total;
}

static void bind_slots(accum *acc)
{
    double *cursor = acc->arena;

    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++) {
        acc->slot[id] = cursor;
        cursor += accum_extent(id, acc->cap);
    }
}

int accum_alloc(accum *acc, size_t cap)
{
    acc->cap   = cap;
    acc->arena = calloc(arena_extent(cap), sizeof *acc->arena);
    if (!acc->arena)
        return -1;

    bind_slots(acc);
    return 0;
}

void accum_free(accum *acc)
{
    free(acc->arena);
    acc->arena = NULL;
}

void accum_zero(accum *acc, size_t len)
{
    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++)
        memset(acc->slot[id], 0, accum_extent(id, len) * sizeof *acc->arena);
}

void accum_add(accum *dst, const accum *src, size_t len)
{
    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++) {
        const double *from = src->slot[id];
        double       *into = dst->slot[id];
        size_t        n    = accum_extent(id, len);

        for (size_t i = 0; i < n; i++)
            into[i] += from[i];
    }
}
