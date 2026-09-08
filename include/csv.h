/* csv.h -- a cmuts output written as comma separated values.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdio.h>

#include "format.h"

typedef struct {
    const char *input_path;  /* the output to convert */
    const char *fasta_path;  /* the FASTA naming the rows, or NULL to number them */
} csv_config;

/* Writes one row per position of every reference to out, with one column for each
 * requested field the input holds. Returns 0, or -1 with a description in error. */
int csv_run(const csv_config *cfg, const fmt_request *requests, size_t n_requests,
            FILE *out, char *error, size_t error_len);
