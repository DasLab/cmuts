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

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

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
    OUT_UNMAPPED,
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

/* One field of the output.
 *
 * A field indexed by reference holds one row per reference, and its dataset is the
 * reference dimension followed by the extents of that row. A field not indexed by one
 * belongs to the run: the unmapped reads align to no reference, so their total is a row
 * with nothing in front of it, which is a dataset of rank zero. */
typedef struct {
    const char *name;
    shape_fn    row;      /* the extents one reference's values occupy */
    bool        per_ref;  /* whether there is one such row per reference */
    out_stored  stored;   /* the type its values are narrowed to */
    out_absent  absent;   /* what a value it was never given means */
    out_combine combine;  /* how two runs' values of it combine */
} out_field;

extern const out_field OUT_FIELDS[OUT_N_FIELDS];

/* The reference dimension, where a field has one, and the extents of its row. */
#define OUT_RANK_MAX (1 + SHAPE_RANK_MAX)

/* Values one field occupies for a reference of len bases, in a file whose longest
 * reference is cap. */
size_t out_values(out_field_id id, size_t len, size_t cap);

/* The widest row of any field, which is what a buffer must hold to take a row of any of
 * them. */
size_t out_widest(size_t cap);

/* Bytes one of a field's values occupies, and the most any field's value occupies. A
 * buffer taking a row of any field is as long as out_widest values of out_widest_bytes
 * each. */
size_t out_stored_bytes(out_field_id id);
size_t out_widest_bytes(void);

/* Writes the dimensions of one field's dataset into dims, which must have room for
 * OUT_RANK_MAX of them, and returns the rank. Rows hold as many values as the field
 * occupies at the longest reference, so every row of a field is the same width whatever
 * its own reference measures. */
int out_dims(out_field_id id, int32_t n_refs, size_t cap, size_t *dims);

/* Dimensions one field's dataset has. */
int out_rank(out_field_id id);
