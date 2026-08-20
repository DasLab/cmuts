/* bamstream.h -- several alignment files read as one.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdint.h>

#include "bam.h"

/* A single pass over several files, delivering their records as if the files had been
 * merged.
 *
 * Records arrive grouped by reference and not by file, so a reference is still finished
 * the moment the stream moves past it and is visited once across every file.
 *
 * A record refers to its reference by an index into its own file's header, so the headers
 * must agree, which this module does not check. */
typedef struct cm_bam_stream cm_bam_stream;

/* Opens every path and reads the first record of each. There must be at least one path.
 * Returns NULL only when out of memory; anything else is reported through
 * cm_bam_stream_error(). */
cm_bam_stream *cm_bam_stream_open(const char *const *paths, size_t n_paths,
                                  const char *fasta_path, int threads);

void cm_bam_stream_close(cm_bam_stream *stream);

/* Returns a description of the stream's failure, or NULL if it has not failed. */
const char *cm_bam_stream_error(const cm_bam_stream *stream);

/* Fills out with the next alignment. Returns a cm_iter_status; on CM_ITER_ERROR the cause is
 * available from cm_bam_stream_error(). The record is invalidated by the next call, as
 * cm_bam_next()'s is. A record carrying BAM_FPAIRED is an error: mates must be merged
 * before alignment, or their overlap counts a molecule twice. */
int cm_bam_stream_next(cm_bam_stream *stream, cm_bam_record *out);

/* Gives the htslib record behind the most recent cm_bam_stream_next(), for callers that
 * must retain a read past the next advance. Copy it with bam_copy1. */
const bam1_t *cm_bam_stream_raw(const cm_bam_stream *stream);

/* Give how far the stream has read, and how far it has to go, in compressed bytes of the
 * alignments only, summed over the files. Span is 0 where the size of any one of them is
 * unknown. */
uint64_t cm_bam_stream_position(const cm_bam_stream *stream);
uint64_t cm_bam_stream_span(const cm_bam_stream *stream);

/* ------------------------------------------------------------------------ */
/* Header queries                                                            */
/* ------------------------------------------------------------------------ */

/* Answered from the first file's header. */
const char *cm_bam_stream_refname(const cm_bam_stream *stream, int32_t tid);
int32_t     cm_bam_stream_nref(const cm_bam_stream *stream);
hts_pos_t   cm_bam_stream_reflen(const cm_bam_stream *stream, int32_t tid);
hts_pos_t   cm_bam_stream_max_reflen(const cm_bam_stream *stream);

/* Give the files themselves, for checks that must reach each header in turn. */
size_t               cm_bam_stream_count(const cm_bam_stream *stream);
const cm_bam_reader *cm_bam_stream_reader(const cm_bam_stream *stream, size_t file);
const char          *cm_bam_stream_path(const cm_bam_stream *stream, size_t file);
