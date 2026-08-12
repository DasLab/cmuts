/* rates.h -- what the evidence gathered at a position comes to.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "accum.h"

/* How much evidence a position must carry before a rate is reported for it.
 *
 * A rate over a fraction of an observation is not reportable: the standard error of a
 * proportion is bounded by a half only where the evidence is a whole observation, and below
 * that it is divided by a fraction and unbounded. */
typedef struct {
    double min_depth;   /* met by reaching it */
} rate_config;

rate_config rate_defaults(void);

/* The mutations at each of len positions over the evidence for them, and the standard error
 * of that rate, written into out.
 *
 * Both are NaN where the evidence falls short, and at the same positions, so one can be
 * read from the other.
 *
 * out is written as the accumulator is read, so it must not be one of the accumulator's own
 * arrays: a position's evidence is read after the position before it has been written. This
 * is what restrict states, and the writer's scratch row is what satisfies it. */
void rate_reactivity(const rate_config *cfg, const accum *acc, size_t len,
                     double *restrict out);
void rate_error(const rate_config *cfg, const accum *acc, size_t len,
                double *restrict out);
