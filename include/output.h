/* output.h -- what an output file holds.
 *
 * Defines the fields an output is made of: their names, the width of a row of each, and
 * the type each is stored as. Everything that reads or writes an output takes the
 * description from here.
 *
 * This is the specification of the file and nothing more. Where a field's values come
 * from is refrow's, the accumulated fields and the written ones not corresponding one to
 * one; how those values are stored is h5layout's. Nothing here refers to HDF5, so the
 * pipeline, the subtraction and everything else working in terms of fields and values
 * compiles without it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

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
    OUT_PAIRWISE_CORRELATION,
    OUT_PAIRWISE_CONDITIONAL,
    OUT_PAIRWISE_COVERAGE,
    OUT_NORM,
    OUT_N_FIELDS,
} out_field_id;

/* What a value the run never wrote means. Counts read zero; rates read NaN, since zero
 * would claim an unread position was measured and unmodified.
 *
 * Columns past a reference's length are NaN regardless of this setting, being outside the
 * reference and not unmeasured within it. */
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
} out_field;

extern const out_field OUT_FIELDS[OUT_N_FIELDS];

/* What an output carries besides its fields.
 *
 * These describe the run and not any reference, so they are written as attributes on the
 * root group: everything reading an output walks its datasets, and an attribute stays out
 * of that. Their values differ from one run to the next, so the table names them and the
 * writer is given what each holds. */
typedef enum {
    OUT_ATTR_PROGRAM,
    OUT_ATTR_VERSION,
    OUT_N_ATTRS,
} out_attr_id;

typedef struct {
    const char *name;
    const char *detail;  /* what it holds, in one sentence */
} out_attribute;

extern const out_attribute OUT_ATTRIBUTES[OUT_N_ATTRS];

/* The reference dimension, where a field has one, and the extents of its row. */
#define OUT_RANK_MAX (1 + SHAPE_RANK_MAX)

/* Gives the values one field occupies for a reference of len bases, in a file whose
 * longest reference is cap. */
size_t out_values(out_field_id id, size_t len, size_t cap);

/* One field of one program's output.
 *
 * A field has the same shape, type and fill wherever it is written, so those stay in
 * OUT_FIELDS. What the numbers mean depends on what produced them -- a rate cmuts-hmm
 * counted is not a rate cmuts-sub took a background off -- so the description belongs
 * here, and falls back to the field's own where the two agree. */
typedef struct {
    out_field_id id;
    const char  *detail;     /* what the numbers are here, in one sentence */
    const char  *condition;  /* what a run needs for it, where it is not written always */
} out_written;

/* The fields one program writes.
 *
 * A program is an interface: it takes the datasets it needs, ignores whatever else its
 * input carries, and writes exactly these. */
typedef struct {
    const out_written *fields;
    size_t             n_fields;
} out_manifest;

/* The fields every output holds, which is what a program reading one takes. */
extern const out_manifest OUT_COMMON;

/* Fills one entry per field with whether the manifest writes it. */
void out_selection(const out_manifest *manifest, bool *wanted);

/* Whether a selection holds a field. */
bool out_wanted(out_field_id id, const bool *wanted);

/* Gives the widest row of any field, which is what a buffer must hold to take a row of
 * any of them. */
size_t out_widest(size_t cap, const bool *wanted);

/* Give the bytes one of a field's values occupies, and the most any field's value
 * occupies. A buffer taking a row of any field is as long as out_widest values of
 * out_widest_bytes each. */
size_t out_stored_bytes(out_field_id id);
size_t out_widest_bytes(void);

/* Writes the dimensions of one field's dataset into dims, which must have room for
 * OUT_RANK_MAX of them, and returns the rank. Rows hold as many values as the field
 * occupies at the longest reference, so every row of a field is the same width whatever
 * its own reference measures. */
int out_dims(out_field_id id, int32_t n_refs, size_t cap, size_t *dims);

/* Gives the dimensions one field's dataset has. */
int out_rank(out_field_id id);

/* Writes the table above as JSON, for generating the documentation of the format from the
 * program that writes it. */
void out_dump_layout(FILE *out, const char *program,
                     const out_manifest *manifest);
