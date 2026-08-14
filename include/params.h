/* params.h -- the pair HMM's rates, as a struct and as a file.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdio.h>

typedef struct {
    double open_deletion;
    double open_insertion;
    double extend_deletion;
    double extend_insertion;
    double modification;
} phmm_params;

phmm_params phmm_defaults(void);

/* Reads an HMM parameter file. Returns 0, or -1 with a reason in error. */
int params_read(const char *path, phmm_params *params, char *error, size_t error_len);

/* Writes the defaults in the form params_read accepts. The command line
 * answers --dump-params with it. */
void params_dump_defaults(FILE *out);
