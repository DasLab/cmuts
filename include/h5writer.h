/* h5writer.h -- per-reference results, written as one row per reference.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "output.h"

/* Writes one dataset per field of the output format: (n_refs, width) for a row field,
 * (n_refs,) for a scalar.
 *
 * A row is identified by its position alone: row i belongs to the reference the header
 * declares i-th, which is the i-th record of the FASTA, the two being required to agree.
 * The names are not written, the FASTA already holding them in that order.
 *
 * Every row exists from the outset, the reference count being known before anything is
 * written, so rows may be written in any order and no dataset is ever extended. Unwritten
 * rows and positions past a reference's length keep the field's fill value, which
 * separates "no reads" and "not part of this reference" from a measured zero.
 *
 * HDF5 is not thread-safe unless built for it, so a writer must be used from one thread
 * only. */
typedef struct h5writer h5writer;

/* Whether a file may be created at this path, reporting through may_replace whether the
 * create may replace a file already there.
 *
 * Called before anything is opened, so that a mistyped path costs nothing and no previous
 * result is at risk while the inputs are still being read. Returns 0, or -1 with a reason
 * in error. */
int h5writer_may_replace(const char *path, bool overwrite, bool *may_replace,
                         char *error, size_t error_len);

/* Fails rather than replacing an existing file unless overwrite is set. */
h5writer *h5writer_create(const char *path, int32_t n_refs, size_t ref_cap,
                          bool overwrite);
void      h5writer_close(h5writer *w);

/* Writes one field's values for a reference of len bases, narrowing the accumulated doubles
 * to the field's stored type. Columns past the reference's own extent are marked as not
 * part of it.
 *
 * Pass ref_cap as len to write a full-width row, at which every field occupies its full
 * width and none has a tail to mark. */
int h5writer_field(h5writer *w, out_field_id id, int32_t tid, size_t len,
                   const double *values);

/* Attaches a run-level total to the file, for counts belonging to no single reference and
 * so having no row of their own. */
int h5writer_count(h5writer *w, const char *name, size_t value);

const char *h5writer_error(const h5writer *w);
