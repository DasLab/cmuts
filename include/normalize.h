/* normalize.h -- reactivity rates divided by a scale taken from the rates themselves.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "output.h"

/* How the scale is taken from the pooled rates. */
typedef enum {
    NORM_UBR,
    NORM_OUTLIER,
} norm_scheme;

typedef struct {
    const char *const *inputs;
    const char *const *outputs;  /* one per input, in the same order */
    size_t             n_files;

    norm_scheme scheme;
    double      min_coverage;  /* the coverage a position needs to join the ubr pool */

    /* NaN where the bound is not applied. */
    double clip_above;

    bool overwrite;
} normalize_config;

/* Divides every input by one scale pooled over all of them, writing each to its own
 * output. program is recorded in each as what produced it. Returns 0, or -1 with a
 * description in error. */
int normalize_run(const normalize_config *cfg, const char *program,
                  const out_manifest *writes, char *error,
                  size_t error_len);
