/* format.h -- the file format every program reads and writes.
 *
 * Defines the fields a file is made of: their names, the width of a row of each, and
 * the type each is stored as. Everything that reads or writes a file takes the
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
    FMT_COVERAGE,
    FMT_REACTIVITY,
    FMT_ERROR,
    FMT_LENGTHS,
    FMT_READS,
    FMT_REJECTED,
    FMT_UNMAPPED,
    FMT_PAIRWISE_CORRELATION,
    FMT_PAIRWISE_CONDITIONAL,
    FMT_PAIRWISE_COVERAGE,
    FMT_NORM,
    FMT_SEQUENCE,
    FMT_N_FIELDS,
} fmt_field_id;

/* The type a field's values are narrowed to in the output file. */
typedef enum {
    FMT_F32,
    FMT_U64,
    FMT_I8,
    FMT_N_STORED,
} fmt_stored;

/* One field of the output. A field indexed by reference holds one row per reference, and
 * its dataset is the reference dimension followed by the extents of that row. */
typedef struct {
    const char *name;
    const char *detail;    /* what the dataset holds, in one sentence */
    shape_fn    row;       /* the extents one reference's values occupy */
    bool        per_ref;   /* whether there is one such row per reference */
    bool        from_ref;  /* whether its values are the reference's and not the reads' */
    fmt_stored  stored;    /* the type its values are narrowed to */
    double      fill;      /* what a position the run never wrote reads as */
} fmt_field;

extern const fmt_field FMT_FIELDS[FMT_N_FIELDS];

/* What an output carries besides its fields, as attributes on the root group. */
typedef enum {
    FMT_ATTR_PROGRAM,
    FMT_ATTR_VERSION,
    FMT_N_ATTRS,
} fmt_attr_id;

typedef struct {
    const char *name;
    const char *detail;  /* what it holds, in one sentence */
} fmt_attribute;

extern const fmt_attribute FMT_ATTRIBUTES[FMT_N_ATTRS];

/* The reference dimension, where a field has one, and the extents of its row. */
#define FMT_RANK_MAX (1 + SHAPE_RANK_MAX)

/* Gives the values one field occupies for a reference of len bases, in a file whose
 * longest reference is cap. */
size_t fmt_values(fmt_field_id id, size_t len, size_t cap);

/* One field of one program's output. */
typedef struct {
    fmt_field_id id;
    /* What this program changes about the field, in one sentence. NULL where the field's
     * own detail says it all. */
    const char  *note;
    const char  *condition;  /* what a run needs for it, where it is not written always */
    /* The fields the run reads to write this one, ending in FMT_N_FIELDS. NULL for a
     * field written without reading any. */
    const fmt_field_id *depends;
    /* Whether a run that cannot write this field is refused. The fields it depends on
     * are then required of every input. */
    bool required;
} fmt_written;

/* Names the fields one entry's depends lists. */
#define FMT_DEPENDS(...) ((const fmt_field_id[]){ __VA_ARGS__, FMT_N_FIELDS })

/* The fields one program writes, and, through each entry's depends, the fields it reads
 * to write them. */
typedef struct {
    const fmt_written *fields;
    size_t             n_fields;
} fmt_manifest;

/* Fills one entry per field with whether the manifest holds it. */
void fmt_selection(const fmt_manifest *manifest, bool *wanted);

/* One field a program asks to read of an input file. */
typedef struct {
    fmt_field_id id;
    bool         required;  /* an input lacking it is refused */
} fmt_request;

/* Fills requests, which must hold FMT_N_FIELDS entries, with every field the manifest's
 * entries depend on, each once, required where a required entry depends on it. Returns
 * how many it filled. */
size_t fmt_requests_of(const fmt_manifest *manifest, fmt_request *requests);

/* Whether a selection holds a field. */
bool fmt_wanted(fmt_field_id id, const bool *wanted);

/* Gives the widest row of any field, which is what a buffer must hold to take a row of
 * any of them. */
size_t fmt_widest(size_t cap, const bool *wanted);

/* Whether this field's values must be written for a reference no read arrived on. */
bool fmt_values_needed(fmt_field_id id);

/* One value of any field, in the type it is stored as. */
typedef union {
    float    f32;
    uint64_t u64;
    int8_t   i8;
} fmt_value;

/* Narrows a field's fill to the type it is stored in. Returns 0, or -1 where that type
 * has no such value: a whole number has no NaN. */
int fmt_fill_value(fmt_field_id id, fmt_value *value);

/* Give the bytes one of a field's values occupies, and the most any field's value
 * occupies. */
size_t fmt_stored_bytes(fmt_field_id id);
size_t fmt_widest_bytes(void);

/* Writes the dimensions of one field's dataset into dims, which must have room for
 * FMT_RANK_MAX of them, and returns the rank. Every row spans the longest reference, so
 * a field's rows are all one width. */
int fmt_dims(fmt_field_id id, int32_t n_refs, size_t cap, size_t *dims);

/* Gives the rank of one field's dataset. */
int fmt_rank(fmt_field_id id);

/* Writes the table above as JSON, for generating the documentation of the format from the
 * program that writes it. */
void fmt_dump_layout(FILE *out, const char *program,
                     const fmt_manifest *manifest);
