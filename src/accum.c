/* accum.c -- accumulator storage, laid out from the field table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "accum.h"

#include <stdlib.h>
#include <string.h>

const accum_field ACCUM_FIELDS[ACCUM_N_FIELDS] = {
    [ACCUM_COVERAGE]  = { "coverage",       ACCUM_PER_BASE },
    [ACCUM_MUTATIONS] = { "mutations",      ACCUM_PER_BASE },
    [ACCUM_READS]     = { "reads",          ACCUM_SCALAR   },
    [ACCUM_FILTERED]  = { "reads_filtered", ACCUM_SCALAR   },
};

/* Values a field occupies when the reference under consideration is len bases
 * long; with len == cap this also gives the field's stride in the arena. */
static size_t field_extent(accum_field_id id, size_t len)
{
    return ACCUM_FIELDS[id].kind == ACCUM_PER_BASE ? len : 1;
}

static size_t arena_extent(size_t len)
{
    size_t total = 0;

    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++)
        total += field_extent(id, len);

    return total;
}

size_t accum_bytes(size_t cap)
{
    return arena_extent(cap) * sizeof(double);
}

static void bind_slots(accum *acc)
{
    double *cursor = acc->arena;

    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++) {
        acc->slot[id] = cursor;
        cursor += field_extent(id, acc->cap);
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
        memset(acc->slot[id], 0, field_extent(id, len) * sizeof *acc->arena);
}

void accum_add(accum *dst, const accum *src, size_t len)
{
    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++) {
        const double *from = src->slot[id];
        double       *into = dst->slot[id];
        size_t        n    = field_extent(id, len);

        for (size_t i = 0; i < n; i++)
            into[i] += from[i];
    }
}
