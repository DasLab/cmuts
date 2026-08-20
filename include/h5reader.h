/* h5reader.h -- reading an output file back, one row at a time.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdint.h>

#include "output.h"

/* Reads the datasets h5writer produces, taking the same layout from output.h. Rows are read
 * one at a time, so memory is bounded by the longest reference and not by the size of the
 * file.
 * Values are returned in the type they are stored as, so a value read here and written
 * back out is never widened on the way.
 * HDF5 is not thread-safe unless built for it, so a reader must be used from one thread
 * only. */
typedef struct h5reader h5reader;

/* Opens a file, taking the fields the manifest names and checking that their shapes agree.
 * A field the manifest requires must be present; one it does not require is skipped where
 * the file lacks it. Returns NULL only when out of memory; every other failure is reported
 * through h5reader_error. */
h5reader *h5reader_open(const char *path, const out_manifest *manifest);
void      h5reader_close(h5reader *r);

/* Give the shape the file was written at: one row per reference, each as wide as the
 * longest reference required. */
int32_t h5reader_refs(const h5reader *r);
size_t  h5reader_capacity(const h5reader *r);

/* Whether the reader opened a field. */
bool h5reader_holds(const h5reader *r, out_field_id id);

/* Reads one field's whole row for a reference. values must hold
 * out_values(id, capacity, capacity) of the field's stored type. */
int h5reader_field(h5reader *r, out_field_id id, int32_t tid, void *values);

/* Reads the whole of a field belonging to the run and not to any one reference. */
int h5reader_total(h5reader *r, out_field_id id, size_t *value);

const char *h5reader_error(const h5reader *r);

/* Writes why the reader failed into error, naming the path it was opened from, and
 * returns -1. */
int h5reader_fail(const h5reader *r, const char *path, char *error, size_t error_len);
