/* accum.h -- per-reference accumulated values.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include "shape.h"

/* The quantities every accumulator carries. This enum and ACCUM_FIELDS are the single
 * source of truth for the layout: allocation, zeroing and merging all derive from them.
 * These are not the fields written out; output.h keeps a table of its own. */
typedef enum {
    ACCUM_COVERAGE,
    ACCUM_SPANNED,
    ACCUM_MUTATIONS,
    ACCUM_LENGTHS,
    ACCUM_READS,
    ACCUM_FILTERED,
    ACCUM_N_FIELDS,
} accum_field_id;

/* What the accumulator needs of a field. */
typedef struct {
    shape_fn shape;
} accum_field;

extern const accum_field ACCUM_FIELDS[ACCUM_N_FIELDS];

/* Gives the values one field occupies for a reference of len bases, in a run whose
 * longest is cap. */
size_t accum_values(accum_field_id id, size_t len, size_t cap);

/* Accumulated values for a single reference.
 *
 * Every accumulator in a run is allocated at the same capacity, so any two are layout
 * compatible and may be merged. Values are held as double and narrowed only on output:
 * one reference can absorb hundreds of thousands of reads, enough for float rounding to
 * become visible in the totals. */
typedef struct {
    double *arena;                  /* owned; one allocation backing every field */
    double *slot[ACCUM_N_FIELDS];   /* views into arena */
    size_t  cap;                    /* capacity, in bases */
} accum;

int  accum_alloc(accum *acc, size_t cap);
void accum_free(accum *acc);

/* Zeroes the first len values of every per-base field, and every scalar. */
void accum_zero(accum *acc, size_t len);

/* dst += src over the first len values of every per-base field, and over every
 * scalar. Both accumulators must have the same capacity. */
void accum_add(accum *dst, const accum *src, size_t len);

/* Give the storage for one field: len values for SHAPE_PER_BASE, one for SHAPE_SCALAR. */
static inline double *accum_data(accum *acc, accum_field_id id)
{
    return acc->slot[id];
}

static inline const double *accum_const_data(const accum *acc, accum_field_id id)
{
    return acc->slot[id];
}
