/* output.h -- what an output file holds.
 *
 * Defines the fields an output is made of: their names, the width of a row of
 * each, and how two files' values of one combine under background subtraction.
 * Everything that reads or writes an output takes the description from here, so
 * it cannot fall out of step between them.
 *
 * How those fields are stored is a separate matter, and is h5layout's. Nothing
 * here names HDF5, so the pipeline, the subtraction and everything else working
 * in terms of fields and values compiles without it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "accum.h"

/* What the output holds, which is not what the accumulator holds. The
 * mutations and the span are the evidence gathered; what a caller wants is the
 * rate they come to and how far it is to be believed, so those are what is
 * written and the two they are made of stay inside. */
typedef enum {
    OUT_COVERAGE,
    OUT_REACTIVITY,
    OUT_ERROR,
    OUT_LENGTHS,
    OUT_READS,
    OUT_REJECTED,
    OUT_N_FIELDS,
} out_field_id;

/* How a field's values from two files are combined under background
 * subtraction.
 *
 * Coverage and the read counts are totals, so they add. Reactivity is a rate,
 * so the background's is subtracted from the treated one. The error of that
 * difference is the two errors added in quadrature, the two runs being
 * independent. */
typedef enum {
    OUT_ADD,
    OUT_SUBTRACT,
    OUT_PROPAGATE,
} out_combine;

/* The meaning of a value the run never wrote.
 *
 * For a count, zero is correct: a position no read reached has a count of zero,
 * and the absence is itself the measurement. A rate is different. A reference
 * with no reads has no rate at all, and zero would mark its every position as
 * measured and found unmodified -- the most confident result the output can
 * hold, for a reference nothing was observed of. Those fields are left NaN,
 * which no arithmetic turns back into a result.
 *
 * Applies to values never written. Columns past a reference's own length are
 * marked NaN regardless of this setting, being outside the reference rather
 * than unmeasured within it, which is what makes coverage a record of a
 * reference's extent as well as of its depth. */
typedef enum {
    OUT_ZERO,
    OUT_NAN,
} out_absent;

/* A count is whole however it was accumulated: the arena is doubles throughout,
 * so that any two accumulators merge, and what is only ever a tally of ones
 * comes back out of it as one. Written as a float it would round above two to
 * the twenty-fourth, which a deeply read reference passes. What is fractional
 * in its own right -- the coverage, the rate, its error -- stays float. */
typedef struct {
    const char    *name;
    accum_field_id shape;    /* the accumulated field whose extent it takes */
    bool           counted;  /* whole, and written as an unsigned */
    out_absent     absent;   /* what a value it was never given means */
    out_combine    combine;  /* how two runs' values of it come together */
} out_field;

extern const out_field OUT_FIELDS[OUT_N_FIELDS];

/* The group holding the counts about reads. The field names above include this
 * prefix, and run totals are named relative to it. */
#define OUT_READS_GROUP "reads"

/* A scalar field is one value per reference; every other kind is a reference
 * by something, and so a row. Both shapes are described through the same
 * arrays, which are sized to the larger. */
#define OUT_RANK_SCALAR 1
#define OUT_RANK_VECTOR 2
#define OUT_RANK_MAX    OUT_RANK_VECTOR

/* Values one field occupies for a reference of len bases, in a file whose
 * longest reference is cap. */
size_t out_extent(out_field_id id, size_t len, size_t cap);

/* The widest row of any field, which is the size a buffer must have to hold a
 * row of any of them. */
size_t out_widest(size_t cap);

int out_rank(out_field_id id);
