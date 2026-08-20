/* output.h -- what an output file holds.
 *
 * Defines the fields an output is made of: their names, the width of a row of each, and
 * the type each is stored as. Everything that reads or writes an output takes the
 * description from here.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "shape.h"

/* The fields an output holds. */
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
    OUT_SEQUENCE,
    OUT_N_FIELDS,
} out_field_id;

/* The type a field's values are narrowed to in the output file. */
typedef enum {
    OUT_F32,
    OUT_U64,
    OUT_I8,
    OUT_N_STORED,
} out_stored;

/* One field of the output. A field indexed by reference holds one row per reference, and
 * its dataset is the reference dimension followed by the extents of that row. */
typedef struct {
    const char *name;
    shape_fn    row;       /* the extents one reference's values occupy */
    bool        per_ref;   /* whether there is one such row per reference */
    bool        from_ref;  /* whether its values are the reference's and not the reads' */
    out_stored  stored;    /* the type its values are narrowed to */
    double      fill;      /* what a position the run never wrote reads as */
} out_field;

extern const out_field OUT_FIELDS[OUT_N_FIELDS];

/* What an output carries besides its fields, as attributes on the root group. */
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

/* Where a program's field comes from. */
typedef enum {
    OUT_IF_PRESENT,  /* read where an input has it, skipped where it does not */
    OUT_REQUIRED,    /* read from every input; one lacking it is refused */
    OUT_MADE,        /* the run makes it, so no input is read for it */
} out_origin;

/* One field of one program's output. */
typedef struct {
    out_field_id id;
    const char  *detail;     /* what the numbers are here, in one sentence */
    const char  *condition;  /* what a run needs for it, where it is not written always */
    out_origin   origin;     /* where the run gets it */
} out_written;

/* The fields one program reads and writes. A program reads exactly these of its inputs,
 * ignores whatever else they carry, and writes exactly these. */
typedef struct {
    const out_written *fields;
    size_t             n_fields;
} out_manifest;

/* Fills one entry per field with whether the manifest holds it. */
void out_selection(const out_manifest *manifest, bool *wanted);

/* Gives where the manifest gets a field, or OUT_IF_PRESENT for one it does not hold. */
out_origin out_origin_of(const out_manifest *manifest, out_field_id id);

/* Whether a selection holds a field. */
bool out_wanted(out_field_id id, const bool *wanted);

/* Gives the widest row of any field, which is what a buffer must hold to take a row of
 * any of them. */
size_t out_widest(size_t cap, const bool *wanted);

/* Whether this field's values must be written for a reference no read arrived on. */
bool out_values_needed(out_field_id id);

/* One value of any field, in the type it is stored as. */
typedef union {
    float    f32;
    uint64_t u64;
    int8_t   i8;
} out_value;

/* Narrows a field's fill to the type it is stored in. Returns 0, or -1 where that type
 * has no such value: a whole number has no NaN. */
int out_fill_value(out_field_id id, out_value *value);

/* Give the bytes one of a field's values occupies, and the most any field's value
 * occupies. */
size_t out_stored_bytes(out_field_id id);
size_t out_widest_bytes(void);

/* Writes the dimensions of one field's dataset into dims, which must have room for
 * OUT_RANK_MAX of them, and returns the rank. Every row spans the longest reference, so
 * a field's rows are all one width. */
int out_dims(out_field_id id, int32_t n_refs, size_t cap, size_t *dims);

/* Gives the rank of one field's dataset. */
int out_rank(out_field_id id);

/* Writes the table above as JSON, for generating the documentation of the format from the
 * program that writes it. */
void out_dump_layout(FILE *out, const char *program,
                     const out_manifest *manifest);
