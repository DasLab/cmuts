/* output.h -- the layout of an output file.
 *
 * Defines the fields an output holds: their names, the width of a row of each,
 * how each is stored, and how two files' values of one combine under background
 * subtraction. The writer, the reader and the subtraction all take the layout
 * from here, so it cannot fall out of step between them.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <hdf5.h>

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

/* A count is whole however it was accumulated: the arena is doubles throughout,
 * so that any two accumulators merge, and what is only ever a tally of ones
 * comes back out of it as one. Written as a float it would round above two to
 * the twenty-fourth, which a deeply read reference passes. What is fractional
 * in its own right -- the coverage, the rate, its error -- stays float. */
typedef struct {
    const char    *name;
    accum_field_id shape;    /* the accumulated field whose extent it takes */
    bool           counted;  /* whole, and written as an unsigned */
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

/* The dataset's shape and that of the chunks it is stored in. The width of a
 * row is the field's own extent at the longest reference, so a field wider than
 * one value per base is sized by the same rule as the rest. Both arrays are
 * OUT_RANK_MAX long, and only the field's own rank is written. */
void out_shape(out_field_id id, int32_t n_refs, size_t cap,
               hsize_t *dims, hsize_t *chunk);

hid_t       out_type(out_field_id id);
const void *out_fill(out_field_id id);

/* A dataspace holding one row of the widest field, which a span of any field's
 * row is then selected from. */
hid_t out_row_space(size_t cap);

/* Selects n values of one reference's row of a field, starting at column from:
 * in the file, and in the memory row the values are moved through. The two are
 * selected together so that they cannot disagree on how many values move.
 *
 * The selections replace whatever was there, so both dataspaces may be kept for
 * as long as the file is open and reselected for each row. */
int out_select_span(hid_t filespace, hid_t memspace, out_field_id id,
                    int32_t tid, size_t from, size_t n);

/* A creation property list with object timestamping turned off.
 *
 * HDF5 stamps every object header with the time it was written, so two runs
 * over the same input would produce files differing in bytes that say nothing
 * about the result. Applied to the file as well as the datasets: an fcpl also
 * carries the root group's creation properties, and that group is stamped like
 * any other. */
hid_t out_untimed_plist(hid_t class_id);

/* A dataset creation property list: chunked, filtered, and filled as the field
 * requires. */
hid_t out_layout_plist(out_field_id id, const hsize_t *chunk, int rank);

/* A dataset access property list whose chunk cache holds several chunks of
 * this shape.
 *
 * The writer and the reader both work through rows in roughly ascending order,
 * so caching a few chunks avoids inflating a chunk again for each row in it. */
hid_t out_access_plist(const hsize_t *chunk, int rank);
