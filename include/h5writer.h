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
 * (n_refs,) for one value per reference, and a scalar for a total belonging to the run.
 *
 * A row is identified by position and not by name: row i belongs to the reference the header
 * declares i-th, which is the i-th record of the FASTA, the two being required to agree.
 * The names are not written, the FASTA already holding them in that order.
 *
 * Every row exists from the outset, the reference count being known before anything is
 * written, so rows may be written in any order and no dataset is ever extended. Unwritten
 * rows and positions past a reference's length keep the field's fill value, which
 * separates "no reads" and "not part of this reference" from a measured zero.
 *
 * HDF5 is not thread-safe unless built for it, so a writer must be used from one thread
 * only -- with the exceptions of h5writer_expect and h5chunk_filter, which touch no HDF5
 * state.
 *
 * A gathering writer collects the fields whose rows hold more than one value into whole
 * chunks and writes each with H5Dwrite_chunk, since filtering rows one at a time on the
 * writing thread otherwise dominates the run once the references number in the thousands.
 * A chunk handed to H5Dwrite_chunk must arrive already filtered, and that work belongs to
 * no particular thread: the writing thread takes each finished chunk with
 * h5writer_take_chunk, any thread runs h5chunk_filter over it, and the writing thread
 * hands it back through h5writer_write_chunk. Every chunk taken must come back before the
 * writer is closed; whatever was never taken, the close filters and writes itself. */
typedef struct h5writer h5writer;

/* One finished chunk of one field, carrying its values and the buffer its filtered bytes
 * go to. It travels between threads without the writer, so it also carries the geometry
 * that filtering needs. */
typedef struct h5chunk h5chunk;

/* Whether a file may be created at this path, reporting through may_replace whether the
 * create may replace a file already there.
 *
 * Called before anything is opened, so that a mistyped path costs nothing and no previous
 * result is at risk while the inputs are still being read. Returns 0, or -1 with a reason
 * in error. */
int h5writer_may_replace(const char *path, bool overwrite, bool *may_replace,
                         char *error, size_t error_len);

/* Creates the file, failing where one is already at the path unless overwrite is set.
 *
 * program names what is writing the file, and is recorded on the root group alongside the
 * version, so that a result read back later can be traced to what produced it.
 *
 * wanted names the optional fields the run writes, one entry per field, and may be NULL
 * for a run writing every field there is. A field left out has no dataset and no buffer
 * sized for it.
 *
 * gather turns on whole-chunk writing for the wide fields; without it every write goes
 * through the row path, with HDF5 filtering on the writing thread. */
h5writer *h5writer_create(const char *path, const char *program, int32_t n_refs,
                          size_t ref_cap, bool overwrite, const bool *wanted,
                          bool gather);

/* The loading thread, as each reference is opened, in ascending order. A chunk cannot be
 * written until every reference opened in it has arrived, and only the loader knows which
 * those are. Does nothing for a writer that does not gather. */
void h5writer_expect(h5writer *w, int32_t tid);

/* The writing thread, once a reference's fields have all been given. Finishes the chunks
 * nothing further can land in, which h5writer_take_chunk then hands out. Does nothing for
 * a writer that does not gather. */
int h5writer_wrote(h5writer *w, int32_t tid);

/* The next finished chunk awaiting a filter, or NULL while none is. The caller owns it
 * until it comes back through h5writer_write_chunk. */
h5chunk *h5writer_take_chunk(h5writer *w);

/* Filters a chunk: shuffle and deflate, as the datasets declare. Touches no HDF5 state
 * and nothing of the writer, so any thread may run it. */
void h5chunk_filter(h5chunk *chunk);

/* Writes a filtered chunk and frees it, whether or not the write succeeds. */
int h5writer_write_chunk(h5writer *w, h5chunk *chunk);

/* Whether the writer holds a field, so a caller may skip deriving values for one that
 * was left out. */
bool h5writer_holds(const h5writer *w, out_field_id id);

/* Writes one reference's block of a field whose row has two extents. The columns and rows
 * past the reference keep the dataset's fill, so only the block itself is written. */
int h5writer_block(h5writer *w, out_field_id id, int32_t tid, size_t len,
                   const double *values);
void      h5writer_close(h5writer *w);

/* Writes one field's values for a reference of len bases, narrowing the accumulated doubles
 * to the field's stored type, and nothing besides: the columns past them are
 * h5writer_pad's.
 *
 * Pass ref_cap as len to write a full-width row. */
int h5writer_field(h5writer *w, out_field_id id, int32_t tid, size_t len,
                   const double *values);

/* Pads one field's row past a reference of len bases, marking those columns as outside it
 * and writing none of its values. Does nothing for a reference as long as the longest,
 * whose row has no room for padding. */
int h5writer_pad(h5writer *w, out_field_id id, int32_t tid, size_t len);

/* Writes one field's whole row for a reference, in the type the field is stored as. values
 * must hold out_values(id, ref_cap, ref_cap) of that type, every column being written. */
int h5writer_row(h5writer *w, out_field_id id, int32_t tid, const void *values);

/* Writes the whole of a field belonging to the run and not to any one reference, and so
 * having no row. */
int h5writer_total(h5writer *w, out_field_id id, size_t value);

/* Writes a field holding one value for the run, narrowed to the type it is stored as.
 * The counterpart for a field holding a count is h5writer_total. */
int h5writer_value(h5writer *w, out_field_id id, double value);


const char *h5writer_error(const h5writer *w);

/* Writes why the writer failed into error, naming the path it was created at. Stands in
 * for a writer that failed without saying why. Returns -1, which is what a caller reporting
 * a failure returns. */
int h5writer_fail(const h5writer *w, const char *path, char *error, size_t error_len);
