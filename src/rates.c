/* rates.c -- the events of each kind at a position over that channel's denominator.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "rates.h"

#include <math.h>

rate_config rate_defaults(void)
{
    return (rate_config){ .min_depth = 1 };
}

/* The accumulated events each channel counts. */
static const accum_field_id EVENTS_OF[RATE_N_CHANNELS] = {
    [RATE_MISMATCHES] = ACCUM_MISMATCHES,
    [RATE_INSERTIONS] = ACCUM_INSERTIONS,
    [RATE_DELETIONS]  = ACCUM_DELETIONS,
};

/* Where each channel's trials sit. A mismatch is tried on every pairing of the
 * position. Each opening is one arm of the decision a pairing makes toward the next
 * base, so its trials are the pairings that decide: the coverage less the reads whose
 * placed span ends there, at the position itself for an insertion and one base 3' of
 * it for a deletion. */
static const struct {
    size_t shift;      /* bases 3' of the position the deciding pairing sits */
    bool   deciding;   /* count only the pairings a decision follows */
} TRIALS_OF[RATE_N_CHANNELS] = {
    [RATE_MISMATCHES] = { 0, false },
    [RATE_INSERTIONS] = { 0, true },
    [RATE_DELETIONS]  = { 1, true },
};

/* Returns a channel's trials at one position, or zero where the deciding pairing falls
 * past the reference. Rounding can push the difference a hair below zero, so it is
 * floored. */
static double depth_at(rate_channel channel, const double *coverage,
                       const double *ends, size_t i, size_t len)
{
    size_t at = i + TRIALS_OF[channel].shift;
    double depth;

    if (at >= len) {
        return 0.0;
    }

    depth = TRIALS_OF[channel].deciding ? coverage[at] - ends[at] : coverage[at];

    return depth > 0.0 ? depth : 0.0;
}

/* Returns whether a position carries enough depth to report on. Some is always required:
 * a depth of zero does not report on a position with none. */
static bool meets_min_depth(double wanted, double depth)
{
    return depth > 0.0 && depth >= wanted;
}

/* Returns whether a position falls within a masked end of the reference. */
static bool masked_at(const rate_config *cfg, size_t i, size_t len)
{
    return i < cfg->nan_5p || len - i <= cfg->nan_3p;
}

/* Returns the events at a position over its depth, held to one. Every event is a share
 * of the depth it is tried against, so the ratio cannot exceed one except by rounding. */
static double ratio_of(double events, double depth)
{
    double rate = depth > 0.0 ? events / depth : 0.0;

    return rate > 1.0 ? 1.0 : rate;
}

/* Returns whether a position reports a rate: it carries enough depth and lies outside
 * the masked ends. */
static bool reported_at(const rate_config *cfg, double depth, size_t i, size_t len)
{
    return meets_min_depth(cfg->min_depth, depth) && !masked_at(cfg, i, len);
}

/* Returns one position's rate, or NaN where none is reported. */
static double rate_at(const rate_config *cfg, double events, double depth,
                      size_t i, size_t len)
{
    if (!reported_at(cfg, depth, i, len)) {
        return (double)NAN;
    }

    return ratio_of(events, depth);
}

/* Returns the binomial standard error of one position's rate, or NaN where none is
 * reported. */
static double error_at(const rate_config *cfg, double events, double depth,
                       size_t i, size_t len)
{
    double rate;

    if (!reported_at(cfg, depth, i, len)) {
        return (double)NAN;
    }

    rate = ratio_of(events, depth);

    return sqrt(rate * (1.0 - rate) / depth);
}

void rate_of(const rate_config *cfg, const accum *acc, size_t len,
             rate_channel channel, double *restrict out)
{
    const double *events   = accum_const_data(acc, EVENTS_OF[channel]);
    const double *coverage = accum_const_data(acc, ACCUM_COVERAGE);
    const double *ends     = accum_const_data(acc, ACCUM_ENDS);

    for (size_t i = 0; i < len; i++) {
        double depth = depth_at(channel, coverage, ends, i, len);

        out[i] = rate_at(cfg, events[i], depth, i, len);
    }
}

void rate_error_of(const rate_config *cfg, const accum *acc, size_t len,
                   rate_channel channel, double *restrict out)
{
    const double *events   = accum_const_data(acc, EVENTS_OF[channel]);
    const double *coverage = accum_const_data(acc, ACCUM_COVERAGE);
    const double *ends     = accum_const_data(acc, ACCUM_ENDS);

    for (size_t i = 0; i < len; i++) {
        double depth = depth_at(channel, coverage, ends, i, len);

        out[i] = error_at(cfg, events[i], depth, i, len);
    }
}
