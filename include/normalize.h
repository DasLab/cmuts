/* normalize.h -- reactivity rates divided by a norm taken from the rates themselves.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "format.h"

/* Where the norm comes from. Every scheme but NORM_VALUE computes it from the pooled
 * rates. */
typedef enum {
    NORM_UBR,
    NORM_OUTLIER,
    NORM_VALUE,
} norm_scheme;

typedef struct {
    const char *const *inputs;
    const char *const *outputs;  /* one per input, in the same order */
    size_t             n_files;

    norm_scheme scheme;
    double      value;         /* the norm itself, under NORM_VALUE */
    double      min_coverage;  /* the coverage a position needs to join the pool */

    bool overwrite;
} normalize_config;

/* Divides every input by one norm, writing each to its own output. program is recorded
 * in each as what produced it. Returns 0, or -1 with a description in error. */
int normalize_run(const normalize_config *cfg, const char *program,
                  const fmt_manifest *writes, char *error,
                  size_t error_len);
