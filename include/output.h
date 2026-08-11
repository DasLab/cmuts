/* output.h -- what an output file holds.
 *
 * Defines the fields an output is made of: their names, the width of a row of each, and
 * how two files' values of one combine under background subtraction. Everything that
 * reads or writes an output takes the description from here.
 *
 * This is the specification of the file and nothing more. Where a field's values come
 * from is refrow's, the accumulated fields and the written ones not corresponding one to
 * one; how those values are stored is h5layout's. Nothing here names HDF5, so the
 * pipeline, the subtraction and everything else working in terms of fields and values
 * compiles without it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include "shape.h"

/* The fields an output holds, which are not those the accumulator holds. The mutations
 * and the span are the evidence gathered; what is written is the rate they come to and
 * its error, the two they are derived from staying internal. */
typedef enum {
    OUT_COVERAGE,
    OUT_REACTIVITY,
    OUT_ERROR,
    OUT_LENGTHS,
    OUT_READS,
    OUT_REJECTED,
    OUT_N_FIELDS,
} out_field_id;

/* How a field's values from two files combine under background subtraction. Coverage and
 * the read counts are totals, so they add. Reactivity is a rate, so the background's is
 * subtracted from the treated one. The error of that difference is the two errors in
 * quadrature, the runs being independent. */
typedef enum {
    OUT_ADD,
    OUT_SUBTRACT,
    OUT_PROPAGATE,
} out_combine;

/* What a value the run never wrote means. Counts read zero; rates read NaN, since zero
 * would claim an unread position was measured and unmodified.
 *
 * Columns past a reference's length are NaN regardless of this setting, being outside the
 * reference rather than unmeasured within it. */
typedef enum {
    OUT_ZERO,
    OUT_NAN,
    OUT_N_ABSENT,
} out_absent;

/* The type a field's values are narrowed to in the output file. */
typedef enum {
    OUT_F32,
    OUT_U64,
    OUT_N_STORED,
} out_stored;

/* One field of the output. */
typedef struct {
    const char *name;
    shape_kind  shape;    /* the kind of row it has */
    out_stored  stored;   /* the type its values are narrowed to */
    out_absent  absent;   /* what a value it was never given means */
    out_combine combine;  /* how two runs' values of it combine */
} out_field;

extern const out_field OUT_FIELDS[OUT_N_FIELDS];

/* The group holding the read counts. The field names above include this prefix, and run
 * totals are named relative to it. */
#define OUT_READS_GROUP "reads"

/* Values one field occupies for a reference of len bases, in a file whose longest
 * reference is cap. */
size_t out_extent(out_field_id id, size_t len, size_t cap);

/* The widest row of any field, which is what a buffer must hold to take a row of any of
 * them. */
size_t out_widest(size_t cap);

/* Dimensions one field's dataset has. */
int out_rank(out_field_id id);
