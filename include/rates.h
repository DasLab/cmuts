/* rates.h -- what the evidence gathered at a position comes to.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "accum.h"

/* How much evidence a position must carry before a rate is reported for it. */
typedef struct {
    double min_depth;   /* met by reaching it */
} rate_config;

rate_config rate_defaults(void);

/* Write into out the mutations at each of len positions over the evidence for them, and
 * the standard error of that rate. Both are NaN where the evidence falls short, and at
 * the same positions. out must not alias the accumulator's arrays, as restrict states. */
void rate_reactivity(const rate_config *cfg, const accum *acc, size_t len,
                     double *restrict out);
void rate_error(const rate_config *cfg, const accum *acc, size_t len,
                double *restrict out);
