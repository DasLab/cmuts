/* accum.h -- per-reference accumulated values.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

/* The quantities every accumulator carries.
 *
 * This enum and ACCUM_FIELDS are the single source of truth for the layout:
 * allocation, zeroing, merging and output all derive from them, and nothing
 * else in the pipeline inspects an accumulator's contents. Adding a quantity
 * means adding one enumerator and one table row. */
typedef enum {
    ACCUM_COVERAGE,
    ACCUM_MUTATIONS,
    ACCUM_READS,
    ACCUM_FILTERED,
    ACCUM_N_FIELDS,
} accum_field_id;

typedef enum {
    ACCUM_PER_BASE,  /* one value per reference base */
    ACCUM_SCALAR,    /* one value per reference */
} accum_kind;

typedef struct {
    const char *name;  /* also the name of the dataset this field is written to */
    accum_kind  kind;
} accum_field;

extern const accum_field ACCUM_FIELDS[ACCUM_N_FIELDS];

/* Accumulated values for a single reference.
 *
 * Every accumulator in a run is allocated at the same capacity, so any two are
 * layout compatible and may be merged. Values are held as double and narrowed
 * only on output: one reference can absorb hundreds of thousands of reads,
 * which is enough for float rounding to become visible in the totals. */
typedef struct {
    double *arena;                  /* owned; one allocation backing every field */
    double *slot[ACCUM_N_FIELDS];   /* views into arena */
    size_t  cap;                    /* capacity, in bases */
} accum;

int  accum_alloc(accum *acc, size_t cap);
void accum_free(accum *acc);

/* Bytes one accumulator of the given capacity occupies. */
size_t accum_bytes(size_t cap);

/* Zeroes the first len values of every per-base field, and every scalar. */
void accum_zero(accum *acc, size_t len);

/* dst += src over the first len values of every per-base field, and over every
 * scalar. Both accumulators must have the same capacity. */
void accum_add(accum *dst, const accum *src, size_t len);

/* Storage for one field: len values for ACCUM_PER_BASE, one for ACCUM_SCALAR. */
static inline double *accum_data(const accum *acc, accum_field_id id)
{
    return acc->slot[id];
}
