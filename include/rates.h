/* rates.h -- what the evidence gathered at a position comes to.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "accum.h"

/* How much evidence a position must carry before it is reported on.
 *
 * A rate over nothing is no rate, and one over almost nothing is a number a
 * caller would have to know to distrust: the standard error of a proportion
 * runs to a half and no further, but only where the evidence is a whole
 * observation. Below that it is divided by a fraction and goes where it likes.
 */
typedef struct {
    double min_depth;   /* met by reaching it */
} rate_config;

rate_config rate_defaults(void);

/* The mutations at each of len positions over the evidence for them, and the
 * standard error of that rate, written into out.
 *
 * Both are NaN where the evidence falls short, which is the one thing a caller
 * cannot mistake for a measurement. Which positions those are is the same for
 * the two, so a caller reading one knows what the other holds.
 *
 * out is written as the accumulator is read, so it must not be one of the
 * accumulator's own arrays: a position's evidence is wanted after the position
 * before it has been written. Somewhere of the caller's own, which is what
 * restrict says and what the writer's scratch row is.
 */
void rate_reactivity(const rate_config *cfg, const accum *acc, size_t len,
                     double *restrict out);
void rate_error(const rate_config *cfg, const accum *acc, size_t len,
                double *restrict out);
