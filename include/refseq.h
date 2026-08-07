/* refseq.h -- reference sequences, delivered in BAM header order.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdint.h>

#include "bamstream.h"
#include "fasta.h"

/* Streams the FASTA in lockstep with the BAM header's reference order, so a
 * reference of any size costs one forward pass and no index.
 *
 * Correctness rests on the two orders agreeing, which every advance checks by
 * comparing the record against the header: its name and length always, and its
 * MD5 wherever the header declares one. A re-ordered FASTA fails immediately. A
 * substituted one fails wherever there is a checksum to catch it -- and where
 * there is none, nothing short of the alignments themselves could tell it from
 * the sequence they were made against.
 *
 * Every file's header is checked, which is also what holds them to each other:
 * two that agree with the FASTA agree with one another. */
typedef struct refseq_source refseq_source;

/* Returns NULL on failure, leaving the reason in why. */
refseq_source *refseq_open(const char *fasta_path, const cm_bam_stream *bams,
                           const char **why);
void           refseq_close(refseq_source *src);

/* The sequence for the given reference index. tid must not go backwards, which
 * holds for any coordinate-sorted file. Returns NULL on mismatch or premature
 * end of file, with the cause available from refseq_error(). */
const cm_fasta_record *refseq_advance(refseq_source *src, int32_t tid);

const char *refseq_error(const refseq_source *src);
