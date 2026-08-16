/* refseq.h -- reference sequences, delivered in BAM header order.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdint.h>

#include "bamstream.h"
#include "fasta.h"

/* Streams the FASTA in step with the BAM header's reference order, so a reference of any
 * size costs one forward pass and no index. Each sequence is checked against every file's
 * header as it is read, which is also what holds the files to each other. */
typedef struct refseq_source refseq_source;

/* Configurable checks. Reference length sizes buffers and so is not optional. */
typedef enum {
    REFSEQ_VERIFY_NAME     = 1 << 0,
    REFSEQ_VERIFY_CHECKSUM = 1 << 1,
} refseq_check;

#define REFSEQ_VERIFY_ALL (REFSEQ_VERIFY_NAME | REFSEQ_VERIFY_CHECKSUM)

/* Returns NULL on failure, leaving the reason in why. */
refseq_source *refseq_open(const char *fasta_path, const cm_bam_stream *bams, int verify,
                           const char **why);
void           refseq_close(refseq_source *src);

/* Gives the sequence for the given reference index. tid must not go backwards, which holds
 * for any coordinate-sorted file. Returns NULL on mismatch or premature end of file, with
 * the cause available from refseq_error(). */
const cm_fasta_record *refseq_advance(refseq_source *src, int32_t tid);

const char *refseq_error(const refseq_source *src);
