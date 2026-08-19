/* score.h -- reactivity measured against a known structure.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

typedef struct {
    const char *input_path;       /* the output being scored */
    const char *fasta_path;       /* names the rows, in the order they were written */
    const char *structures_path;  /* dot bracket records, matched by name */
    const char *reference;        /* one reference to score; NULL scores every one */
    const char *bases;            /* the bases a reagent reports on; NULL takes all */
    double      min_coverage;     /* reads a position needs before it is scored */
} score_config;

/* Writes one row per reference and a last row averaging them. Returns 0, or -1 with a
 * description in error. */
int score_run(const score_config *cfg, char *error, size_t error_len);
