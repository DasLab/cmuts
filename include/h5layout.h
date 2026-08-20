/* h5layout.h -- how the fields of an output are stored in HDF5.
 *
 * The other half of output.h: that defines the fields, this defines what they become on
 * disk -- their datatypes, their fill values, the shape of a dataset and of the chunks it
 * is written in, and the dataspaces a row moves through. The writer and the reader share
 * it, so a file written here is described the same way when it is read back.
 *
 * This is where HDF5 enters the library. h5writer.c and h5reader.c are the only sources
 * that include it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdint.h>

#include <hdf5.h>

#include "output.h"

/* The dataset's shape and that of the chunks it is stored in. A row holds as many values
 * as the field occupies at the longest reference, so a field wider than one value per base
 * is sized by the same rule as the rest. Both arrays are OUT_RANK_MAX long, of which only
 * the field's own rank is written. */
void h5layout_shape(out_field_id id, int32_t n_refs, size_t cap,
                    hsize_t *dims, hsize_t *chunk);

hid_t h5layout_type(out_field_id id);
hid_t h5layout_memory_type(out_field_id id);

/* Gives a dataspace holding one row of the widest field, from which a span of any field's
 * row is then selected. */
hid_t h5layout_row_space(size_t cap, const bool *wanted);

/* Selects n values of one reference's row of a field, starting at column from, in the file
 * and in the memory row the values move through. The two are selected together so that
 * they cannot disagree on how many values move.
 *
 * The selections replace whatever was there, so both dataspaces may be kept for as long as
 * the file is open and reselected for each row. */
int h5layout_select_span(hid_t filespace, hid_t memspace, out_field_id id,
                         int32_t tid, size_t from, size_t n);

/* Selects one reference's whole block of a field whose row has two extents, in the file
 * and in memory together. The block is square and starts at the origin of the row, so a
 * reference shorter than the longest leaves the rest of its row unwritten, which the
 * dataset's fill already holds. */
int h5layout_select_block(hid_t filespace, hid_t memspace, out_field_id id,
                          int32_t tid, size_t len);

/* Gives a creation property list with object timestamping turned off.
 *
 * HDF5 stamps every object header with the time it was written, so two runs over the same
 * input would produce files differing in bytes that carry no information about the result.
 * Applied to the file as well as the datasets: an fcpl also carries the root group's
 * creation properties, and that group is stamped like any other. */
hid_t h5layout_untimed_plist(hid_t class_id);

/* Gives a dataset creation property list: chunked, filtered, and filled as the field
 * requires. */
hid_t h5layout_creation_plist(out_field_id id, const hsize_t *chunk, int rank);

/* Gives a dataset access property list whose chunk cache holds several chunks of this
 * shape. The writer and the reader both work through rows in roughly ascending order, so
 * caching a few chunks avoids inflating a chunk again for each row in it. */
hid_t h5layout_access_plist(out_field_id id, const hsize_t *chunk, int rank);
