/* refseq.h -- reference sequences, delivered in BAM header order.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdint.h>

#include "bam.h"
#include "fasta.h"

/* Streams the FASTA in lockstep with the BAM header's reference order, so a
 * reference of any size costs one forward pass and no index.
 *
 * Correctness rests on the two orders agreeing, which every advance checks by
 * comparing the record's name and length against the header. A re-ordered or
 * substituted FASTA fails immediately rather than silently scoring reads
 * against the wrong sequence. */
typedef struct refseq_source refseq_source;

refseq_source *refseq_open(const char *fasta_path, const cm_bam_reader *reader);
void           refseq_close(refseq_source *src);

/* The sequence for the given reference index. tid must not go backwards, which
 * holds for any coordinate-sorted file. Returns NULL on mismatch or premature
 * end of file, with the cause available from refseq_error(). */
const cm_fasta_record *refseq_advance(refseq_source *src, int32_t tid);

const char *refseq_error(const refseq_source *src);
