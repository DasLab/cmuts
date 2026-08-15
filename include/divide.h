/* divide.h -- reactivity rates normalized against a denatured control.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

typedef struct {
    const char *rates_path;
    const char *control_path;
    const char *output_path;
    bool        overwrite;
} divide_config;

/* Divides the rates by the control, field by field, into a file of the same layout.
 * Returns 0, or -1 with a description in error. */
int divide_run(const divide_config *cfg, char *error, size_t error_len);
