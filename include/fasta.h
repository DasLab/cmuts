/* fasta.h -- sequential iteration over the records of a FASTA file.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include "iter.h"

/* A single FASTA record.
 *
 * Every pointer borrows memory owned by the reader that produced the record,
 * and is invalidated by the next cm_fasta_next() or cm_fasta_close() call on
 * that reader. Copy anything that must outlive the current iteration step. */
typedef struct {
    const char *name;     /* record ID, up to the first whitespace */
    const char *comment;  /* rest of the header line; NULL when absent */
    const char *seq;      /* NUL-terminated, with line breaks removed */
    size_t      len;      /* length of seq */
} cm_fasta_record;

typedef struct cm_fasta_reader cm_fasta_reader;

/* Opens path for reading; plain, gzip and bgzf inputs are all accepted. No
 * .fai index is required, consulted or created. Returns NULL on failure,
 * leaving the reason in why. */
cm_fasta_reader *cm_fasta_open(const char *path, const char **why);

/* Fills out with the next record. Returns a cm_iter_status; on CM_ITER_ERROR
 * the cause is available from cm_fasta_error(). */
int cm_fasta_next(cm_fasta_reader *reader, cm_fasta_record *out);

void cm_fasta_close(cm_fasta_reader *reader);

/* Returns a description of the reader's failure, or NULL if it has not failed. */
const char *cm_fasta_error(const cm_fasta_reader *reader);
