/* rates.h -- what the evidence gathered at a position comes to.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "accum.h"

/* Which positions get a rate: those with enough evidence and outside the masked ends. */
typedef struct {
    double min_depth;   /* evidence a position must carry, met by reaching it */
    size_t nan_5p;      /* bases at the 5' end written as NaN */
    size_t nan_3p;      /* bases at the 3' end written as NaN */
} rate_config;

rate_config rate_defaults(void);

/* Write into out the mutations at each of len positions over the evidence for them, and
 * the standard error of that rate. Both are NaN at the same positions: where the evidence
 * falls short, and within the masked ends. out must not alias the accumulator's arrays,
 * as restrict states. */
void rate_reactivity(const rate_config *cfg, const accum *acc, size_t len,
                     double *restrict out);
void rate_error(const rate_config *cfg, const accum *acc, size_t len,
                double *restrict out);
