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

/* Returns a channel's trials at one position. Only a deletion is tried at the positions
 * its events fall on, since a deleted base is read by no read. */
static double depth_at(rate_channel channel, double coverage, double events)
{
    return channel == RATE_DELETIONS ? coverage + events : coverage;
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

    for (size_t i = 0; i < len; i++) {
        double depth = depth_at(channel, coverage[i], events[i]);

        out[i] = rate_at(cfg, events[i], depth, i, len);
    }
}

void rate_error_of(const rate_config *cfg, const accum *acc, size_t len,
                   rate_channel channel, double *restrict out)
{
    const double *events   = accum_const_data(acc, EVENTS_OF[channel]);
    const double *coverage = accum_const_data(acc, ACCUM_COVERAGE);

    for (size_t i = 0; i < len; i++) {
        double depth = depth_at(channel, coverage[i], events[i]);

        out[i] = error_at(cfg, events[i], depth, i, len);
    }
}
