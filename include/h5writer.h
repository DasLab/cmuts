/* h5writer.h -- per-reference results, written as one row per reference.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdint.h>

#include "accum.h"

/* Writes one dataset per accumulator field: (n_refs, ref_cap) for per-base
 * fields, (n_refs,) for scalars, alongside the reference names.
 *
 * Every row exists from the outset, since the reference count is known from
 * the BAM header, so rows may be written in any order and no dataset ever has
 * to be extended. Unwritten rows and positions past a reference's own length
 * keep a fill value of NaN, which separates "no reads" and "not part of this
 * reference" from a measured zero.
 *
 * HDF5 is not thread-safe unless built for it, so a writer must be used from
 * one thread only. */
typedef struct h5writer h5writer;

h5writer *h5writer_create(const char *path, int32_t n_refs, size_t ref_cap);
void      h5writer_close(h5writer *w);

/* Records the name of every reference, in header order. */
int h5writer_names(h5writer *w, const char *const *names, int32_t n_refs);

/* Writes one reference's row, narrowing the accumulated doubles to float.
 * Values past len are left at the fill value. */
int h5writer_row(h5writer *w, int32_t tid, size_t len, const accum *acc);

/* Run-level totals, attached to the file root. Unmapped reads belong to no
 * reference and so have no row of their own. */
int h5writer_counts(h5writer *w, size_t reads_total, size_t reads_unmapped);

const char *h5writer_error(const h5writer *w);
