/* score.h -- reactivity measured against a known structure.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

/* The bases a reagent modifies, as the set of them. T and U name one base, so a structure
 * written as RNA and a reference written as DNA give the same bit. */
typedef enum {
    SCORE_BASE_A = 1 << 0,
    SCORE_BASE_C = 1 << 1,
    SCORE_BASE_G = 1 << 2,
    SCORE_BASE_U = 1 << 3,
} score_base;

#define SCORE_BASES_ALL (SCORE_BASE_A | SCORE_BASE_C | SCORE_BASE_G | SCORE_BASE_U)

typedef struct {
    const char *input_path;       /* the output being scored */
    const char *fasta_path;       /* names the rows, in the order they were written */
    const char *structures_path;  /* dot bracket records, matched by name */
    int         bases;            /* score_base bits; a base outside them is not scored */
    double      min_coverage;     /* reads a position needs before it is scored */
} score_config;

/* Writes one row per reference and a last row averaging them. Returns 0, or -1 with a
 * description in error. */
int score_run(const score_config *cfg, char *error, size_t error_len);
