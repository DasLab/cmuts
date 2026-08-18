/* accum.c -- accumulator storage, laid out from the field table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "accum.h"

#include <stdlib.h>
#include <string.h>

/* Coverage and span differ over what a read did not read: a deleted position is spanned
 * but not covered, since no base was read there.
 *
 * The span is the evidence the mutations are taken against. A pairing and a deletion
 * contribute to it whole. An insertion contributes only as far as --insertion-weight
 * makes it a modification, entering the span weighted exactly as it enters the
 * mutations: weighted into one and not the other, it would count as evidence against a
 * modification, and not the absence of evidence. */
const accum_field ACCUM_FIELDS[ACCUM_N_FIELDS] = {
    [ACCUM_COVERAGE]  = { shape_per_base },
    [ACCUM_SPANNED]   = { shape_per_base },
    [ACCUM_MUTATIONS] = { shape_per_base },
    [ACCUM_LENGTHS]   = { shape_per_length },
    [ACCUM_READS]     = { shape_none },
    [ACCUM_FILTERED]  = { shape_none },
};

/* Values one field occupies. With len == cap this is also its stride in the arena. */
size_t accum_values(accum_field_id id, size_t len, size_t cap)
{
    return shape_values(ACCUM_FIELDS[id].shape, len, cap);
}

static size_t arena_values(size_t cap)
{
    size_t total = 0;

    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++) {
        total += accum_values(id, cap, cap);
    }

    return total;
}

static void bind_slots(accum *acc)
{
    double *cursor = acc->arena;

    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++) {
        acc->slot[id] = cursor;
        cursor += accum_values(id, acc->cap, acc->cap);
    }
}

int accum_alloc(accum *acc, size_t cap)
{
    acc->cap   = cap;
    acc->arena = calloc(arena_values(cap), sizeof *acc->arena);
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
               accum_values(id, len, acc->cap) * sizeof *acc->arena);
    }
}

void accum_add(accum *dst, const accum *src, size_t len)
{
    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++) {
        const double *from = src->slot[id];
        double       *into = dst->slot[id];
        size_t        n    = accum_values(id, len, dst->cap);

        for (size_t i = 0; i < n; i++) {
            into[i] += from[i];
        }
    }
}
