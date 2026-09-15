/* reads.h -- examining sequencing reads, and passing them to an aligner.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

typedef enum {
    READS_FASTQ,
    READS_BAM,
    READS_NEITHER,  /* some other format, or a file that could not be opened */
} reads_format;

/* Returns the format of the reads in path. A compressed file is decompressed before its
 * contents are examined, so a gzipped FASTQ gives READS_FASTQ. */
reads_format reads_format_of(const char *path);

/* Examines a BAM file.
 *
 * Sets aligned to true if the header names references, which makes the file a set of
 * alignments. Sets paired to true if the first record is one of a pair. A file that holds
 * no records is not paired. Returns 0, or -1 with a description in error. */
int reads_inspect_bam(const char *path, bool *aligned, bool *paired,
                      char *error, size_t error_len);

/* Writes every record of a BAM file to fd as FASTQ, then closes fd.
 *
 * The base qualities are copied unchanged, because cmuts hmm reads each one as the
 * probability that the base is wrong. A record that holds no base qualities is therefore
 * refused. Returns 0, or -1 with a description in error. */
int reads_write_fastq(const char *path, int fd, char *error, size_t error_len);
