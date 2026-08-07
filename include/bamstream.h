/* bamstream.h -- several alignment files read as one.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdint.h>

#include "bam.h"

/* A single pass over several files, delivering their records as if the files
 * had been merged.
 *
 * Records arrive grouped by reference rather than by file, so a reference is
 * still finished the moment the stream moves past it. Running the files one
 * after another would visit each reference once per file, and the result
 * written for it would be overwritten by the next file's.
 *
 * A record names its reference by an index into its own file's header, so the
 * headers must agree. Nothing here checks that: refseq compares each of them
 * against the FASTA, which says the same thing wherever it matters. */
typedef struct cm_bam_stream cm_bam_stream;

/* Opens every path and reads the first record of each, which is why the
 * reference and the thread pool are settled here. There must be at least one.
 * Returns NULL only when out of memory; anything else is reported through
 * cm_bam_stream_error(). */
cm_bam_stream *cm_bam_stream_open(const char *const *paths, size_t n_paths,
                                  const char *fasta_path, int threads);

void cm_bam_stream_close(cm_bam_stream *stream);

/* Description of the stream's failure, or NULL if it has not failed. */
const char *cm_bam_stream_error(const cm_bam_stream *stream);

/* Fills out with the next alignment. Returns a cm_iter_status; on
 * CM_ITER_ERROR the cause is available from cm_bam_stream_error(). The record
 * is invalidated by the next call, as cm_bam_next()'s is. */
int cm_bam_stream_next(cm_bam_stream *stream, cm_bam_record *out);

/* The htslib record behind the most recent cm_bam_stream_next(), for callers
 * that need to retain a read past the next advance. Copy it with bam_copy1. */
const bam1_t *cm_bam_stream_raw(const cm_bam_stream *stream);

/* How far the stream has got, and how far it has to go, in compressed bytes of
 * the alignments alone, summed over the files. Span is 0 where the size of any
 * of them is not known. */
uint64_t cm_bam_stream_position(const cm_bam_stream *stream);
uint64_t cm_bam_stream_span(const cm_bam_stream *stream);

/* ------------------------------------------------------------------------ */
/* Header queries                                                            */
/* ------------------------------------------------------------------------ */

/* Answered from the first file, the others being held to it by the checks the
 * FASTA makes on each. */
const char *cm_bam_stream_refname(const cm_bam_stream *stream, int32_t tid);
int32_t     cm_bam_stream_nref(const cm_bam_stream *stream);
hts_pos_t   cm_bam_stream_max_reflen(const cm_bam_stream *stream);

/* The files themselves, for the checks that must reach each header in turn. */
size_t               cm_bam_stream_count(const cm_bam_stream *stream);
const cm_bam_reader *cm_bam_stream_reader(const cm_bam_stream *stream, size_t file);
const char          *cm_bam_stream_path(const cm_bam_stream *stream, size_t file);
