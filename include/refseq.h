/* refseq.h -- reference sequences, delivered in BAM header order.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdint.h>

#include "bamstream.h"
#include "fasta.h"

/* Streams the FASTA in lockstep with the BAM header's reference order, so that a reference
 * of any size costs one forward pass and no index.
 *
 * Correctness rests on the two orders agreeing, which every advance checks by comparing the
 * record against the header: its name and length, and its MD5 wherever the header declares
 * one. A re-ordered FASTA fails immediately, and a substituted one fails wherever there is a
 * checksum to catch it.
 *
 * Every file's header is checked, which is also what holds the files to each other: two that
 * agree with the FASTA agree with one another. */
typedef struct refseq_source refseq_source;

/* The comparisons an advance makes. Dropping one leaves the FASTA paired with the header by
 * position alone, so a name no longer checked is a re-ordering no longer caught unless the
 * lengths or the checksums differ. */
typedef enum {
    REFSEQ_VERIFY_NAME     = 1 << 0,
    REFSEQ_VERIFY_LENGTH   = 1 << 1,
    REFSEQ_VERIFY_CHECKSUM = 1 << 2,
} refseq_check;

#define REFSEQ_VERIFY_ALL \
    (REFSEQ_VERIFY_NAME | REFSEQ_VERIFY_LENGTH | REFSEQ_VERIFY_CHECKSUM)

/* Returns NULL on failure, leaving the reason in why. */
refseq_source *refseq_open(const char *fasta_path, const cm_bam_stream *bams, int verify,
                           const char **why);
void           refseq_close(refseq_source *src);

/* Gives the sequence for the given reference index. tid must not go backwards, which
 * holds for any
 * coordinate-sorted file. Returns NULL on mismatch or premature end of file, with the cause
 * available from refseq_error(). */
const cm_fasta_record *refseq_advance(refseq_source *src, int32_t tid);

const char *refseq_error(const refseq_source *src);
