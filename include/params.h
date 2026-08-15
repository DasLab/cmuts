/* params.h -- the pair HMM's rates, as a struct and as a file.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdio.h>

/* The rates the pair HMM runs on. What a transition rate leaves over goes to the step back
 * to a match, so the two opening rates sum below one and each extending rate is below
 * one. */
typedef struct {
    double open_deletion;     /* a match steps to a deletion */
    double open_insertion;    /* a match steps to an insertion */
    double extend_deletion;   /* a deletion runs on */
    double extend_insertion;  /* an insertion runs on */
    double modification;      /* a reference base was modified before it was read */
} phmm_params;

/* The default rates a run uses. Overridden by fields present in the file passed to
 * --params. */
phmm_params phmm_defaults(void);

/* Reads an HMM parameter file. Returns 0, or -1 with a reason in error. */
int params_read(const char *path, phmm_params *params, char *error, size_t error_len);

/* Writes the defaults in the form params_read accepts. The command line
 * answers --dump-params with it. */
void params_dump_defaults(FILE *out);
