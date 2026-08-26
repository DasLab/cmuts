/* rates.h -- each event kind's rate at a position, over its own denominator.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "accum.h"

/* Which positions get a rate: those whose denominator is deep enough and which lie
 * outside the masked ends. */
typedef struct {
    double min_depth;   /* depth a channel's denominator must reach */
    size_t nan_5p;      /* bases at the 5' end written as NaN */
    size_t nan_3p;      /* bases at the 3' end written as NaN */
} rate_config;

rate_config rate_defaults(void);

/* The event kinds a rate is computed for. A mismatch is tried wherever a read base
 * pairs, so its rate is over the coverage. Each opening is one arm of the decision a
 * pairing makes toward the next base, and the pairing a read's placed span ends on
 * makes no decision, so an opening rate is over the coverage less those ends: at the
 * base itself for an insertion, and one base 3' of it for a deletion. */
typedef enum {
    RATE_MISMATCHES,
    RATE_INSERTIONS,
    RATE_DELETIONS,
    RATE_N_CHANNELS,
} rate_channel;

/* Write into out one channel's events at each of len positions over its denominator, and
 * the binomial standard error of that rate. Both are NaN at the same positions: where
 * the denominator falls short, and within the masked ends. out must not alias the
 * accumulator's arrays, as restrict states. */
void rate_of(const rate_config *cfg, const accum *acc, size_t len,
             rate_channel channel, double *restrict out);
void rate_error_of(const rate_config *cfg, const accum *acc, size_t len,
                   rate_channel channel, double *restrict out);
