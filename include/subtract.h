/* subtract.h -- one output read against another taken as its background.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

typedef struct {
    const char *treated_path;
    const char *untreated_path;
    /* NULL where the run has no control, leaving the reactivity a difference of rates
     * rather than a ratio over one. */
    const char *denatured_path;
    const char *output_path;
    bool        overwrite;
    bool        clip;  /* whether a negative reactivity is raised to zero */
} subtract_config;

/* Combines the inputs field by field into a file of the same layout, under the rule each
 * field follows. Returns 0, or -1 with a description in error. */
int subtract_run(const subtract_config *cfg, char *error, size_t error_len);

/* Writes the rule each field follows as JSON, with and without a control, for
 * generating the documentation of the subtraction from the program doing it. */
void subtract_dump_rules(FILE *out);
