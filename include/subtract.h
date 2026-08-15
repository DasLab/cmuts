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
    const char *output_path;
    bool        overwrite;
    bool        clip;  /* whether a negative reactivity is raised to zero */
} subtract_config;

/* Takes the background off the treated run, field by field, into a file of the same
 * layout. program is recorded in the output as what produced it. Returns 0, or -1 with a
 * description in error. */
int subtract_run(const subtract_config *cfg, const char *program, char *error,
                 size_t error_len);
