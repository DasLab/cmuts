/* accum.c -- accumulator storage, laid out from the field table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "accum.h"

#include <stdlib.h>
#include <string.h>

/* Coverage and span differ over what a read did not read: a deleted position is spanned
 * but not covered, and a poorly read one is spanned whole but covered in proportion to
 * the confidence in its base.
 *
 * The span is the evidence the mutations are taken against. A pairing and a deletion
 * contribute to it whole. An insertion contributes only as far as --insertion-weight
 * makes it a modification, entering the span weighted exactly as it enters the
 * mutations: weighted into one and not the other, it would count as evidence against a
 * modification rather than no evidence either way. */
const accum_field ACCUM_FIELDS[ACCUM_N_FIELDS] = {
    [ACCUM_COVERAGE]  = { ACCUM_PER_BASE },
    [ACCUM_SPANNED]   = { ACCUM_PER_BASE },
    [ACCUM_MUTATIONS] = { ACCUM_PER_BASE },
    [ACCUM_LENGTHS]   = { ACCUM_PER_LENGTH },
    [ACCUM_READS]     = { ACCUM_SCALAR },
    [ACCUM_FILTERED]  = { ACCUM_SCALAR },
};

/* Values one field occupies. With len == cap this is also its stride in the arena. */
size_t accum_extent(accum_field_id id, size_t len, size_t cap)
{
    switch (ACCUM_FIELDS[id].kind) {
        case ACCUM_PER_BASE:   return len;
        case ACCUM_PER_LENGTH: return ACCUM_LENGTH_BINS(cap);
        default:               return 1;
    }
}

static size_t arena_extent(size_t cap)
{
    size_t total = 0;

    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++) {
        total += accum_extent(id, cap, cap);
    }

    return total;
}

static void bind_slots(accum *acc)
{
    double *cursor = acc->arena;

    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++) {
        acc->slot[id] = cursor;
        cursor += accum_extent(id, acc->cap, acc->cap);
    }
}

int accum_alloc(accum *acc, size_t cap)
{
    acc->cap   = cap;
    acc->arena = calloc(arena_extent(cap), sizeof *acc->arena);
    if (!acc->arena) {
        return -1;
    }

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
    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++) {
        memset(acc->slot[id], 0,
               accum_extent(id, len, acc->cap) * sizeof *acc->arena);
    }
}

void accum_add(accum *dst, const accum *src, size_t len)
{
    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++) {
        const double *from = src->slot[id];
        double       *into = dst->slot[id];
        size_t        n    = accum_extent(id, len, dst->cap);

        for (size_t i = 0; i < n; i++) {
            into[i] += from[i];
        }
    }
}
