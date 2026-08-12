/* subtract.h -- one output read against another taken as its background.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

typedef struct {
    const char *treated_path;
    const char *untreated_path;
    /* A sample denatured before it was treated, measuring what the reagent does where
     * structure is absent. NULL where the run has none, which leaves the reactivity a
     * difference of rates rather than a ratio over that measurement. */
    const char *denatured_path;
    const char *output_path;
    bool        overwrite;
    bool        clip;  /* whether a negative reactivity is raised to zero */
} subtract_config;

/* Combines the inputs field by field into a file of the same layout, under the rule each
 * field follows. Returns 0, or -1 with a description in error. */
int subtract_run(const subtract_config *cfg, char *error, size_t error_len);
