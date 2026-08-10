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
    ACCUM_SPANNED,
    ACCUM_MUTATIONS,
    ACCUM_LENGTHS,
    ACCUM_READS,
    ACCUM_FILTERED,
    ACCUM_N_FIELDS,
} accum_field_id;

typedef enum {
    ACCUM_PER_BASE,    /* one value per reference base */
    ACCUM_PER_LENGTH,  /* one value per read length, plus an overflow bin */
    ACCUM_SCALAR,      /* one value per reference */
} accum_kind;

/* Bins a read-length histogram covers, given the longest reference in the run:
 * one for every length from 0 to twice it. The range reaches past a reference
 * because a read carrying insertions or soft-clipped ends is longer than the
 * one it aligns to, which on the libraries measured is where most of them fall.
 *
 * Every row is this wide, whatever its own reference measures. A read length is
 * not a position in a reference, so a column that a short reference has no
 * reads for is a count of zero and not padding; sizing rows individually would
 * leave the same column meaning different things in different rows, and reads
 * uncounted where the array had room for them all along.
 *
 * A read longer than the range is counted in no bin. How many there were is the
 * reads total less the row's own sum. */
#define ACCUM_LENGTH_BINS(cap) (2 * (cap) + 1)

typedef struct {
    const char *name;  /* also the name of the dataset this field is written to */
    accum_kind  kind;
} accum_field;

extern const accum_field ACCUM_FIELDS[ACCUM_N_FIELDS];

/* Values one field occupies for a reference of len bases, in a run whose
 * longest is cap. Allocation, zeroing, merging, the output row and the
 * dataset's own width all derive from this, so a new kind of field is described
 * here and nowhere else. */
size_t accum_extent(accum_field_id id, size_t len, size_t cap);

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

/* Zeroes the first len values of every per-base field, and every scalar. */
void accum_zero(accum *acc, size_t len);

/* dst += src over the first len values of every per-base field, and over every
 * scalar. Both accumulators must have the same capacity. */
void accum_add(accum *dst, const accum *src, size_t len);

/* Storage for one field: len values for ACCUM_PER_BASE, one for ACCUM_SCALAR.
 *
 * Two forms rather than one taking a const accumulator and handing back a
 * writable pointer, which would let anything holding a finished reference
 * quietly alter it. */
static inline double *accum_data(accum *acc, accum_field_id id)
{
    return acc->slot[id];
}

static inline const double *accum_const_data(const accum *acc, accum_field_id id)
{
    return acc->slot[id];
}
