/* rates.c -- the mutations at a position against the evidence for them.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "rates.h"

#include <math.h>

rate_config rate_defaults(void)
{
    return (rate_config){ .min_depth = 1 };
}

/* Returns whether a position carries enough evidence to report on. Some is always
 * required: a depth of zero does not report on a position with none. */
static bool meets_min_depth(double wanted, double evidence)
{
    return evidence > 0.0 && evidence >= wanted;
}

/* Returns whether a position falls within a masked end of the reference. */
static bool masked_at(const rate_config *cfg, size_t i, size_t len)
{
    return i < cfg->nan_5p || len - i <= cfg->nan_3p;
}

/* Returns the mutations at a position over the evidence for them, held to one. Every
 * weight is a share of an event and a deletion or insertion spans what it contributes,
 * so the ratio cannot exceed one except by rounding. */
static double rate_of(double mutations, double evidence)
{
    double rate = evidence > 0.0 ? mutations / evidence : 0.0;

    return rate > 1.0 ? 1.0 : rate;
}

/* Returns whether a position reports a rate: it carries enough evidence and lies
 * outside the masked ends. */
static bool reported_at(const rate_config *cfg, double evidence, size_t i, size_t len)
{
    return meets_min_depth(cfg->min_depth, evidence) && !masked_at(cfg, i, len);
}

/* Returns one position's rate, or NaN where none is reported. */
static double reactivity_at(const rate_config *cfg, double mutations, double evidence,
                            size_t i, size_t len)
{
    if (!reported_at(cfg, evidence, i, len)) {
        return (double)NAN;
    }

    return rate_of(mutations, evidence);
}

/* Returns what the reads left undecided at a position: the posterior variance of its
 * count, from the two accumulated moments. Each read's event is Bernoulli with its
 * posterior chance m, contributing m(1 - m); certain calls contribute nothing, so hard
 * counts give zero. Rounding on weighted events can push the difference below zero,
 * which is held there. */
static double ambiguity_of(double mutations, double squared)
{
    double ambiguity = mutations - squared;

    return ambiguity > 0.0 ? ambiguity : 0.0;
}

/* Returns the standard error of one position's rate, or NaN where none is reported.
 * Its variance is the molecular sampling of the rate over the evidence, plus the
 * ambiguity of the calls behind it; with every call certain the second term vanishes
 * and the plain binomial error remains. */
static double error_at(const rate_config *cfg, double mutations, double squared,
                       double evidence, size_t i, size_t len)
{
    double rate;

    if (!reported_at(cfg, evidence, i, len)) {
        return (double)NAN;
    }

    rate = rate_of(mutations, evidence);

    return sqrt(rate * (1.0 - rate) / evidence
              + ambiguity_of(mutations, squared) / (evidence * evidence));
}

void rate_reactivity(const rate_config *cfg, const accum *acc, size_t len,
                     double *restrict out)
{
    const double *evidence  = accum_const_data(acc, ACCUM_SPANNED);
    const double *mutations = accum_const_data(acc, ACCUM_MUTATIONS);

    for (size_t i = 0; i < len; i++) {
        out[i] = reactivity_at(cfg, mutations[i], evidence[i], i, len);
    }
}

void rate_error(const rate_config *cfg, const accum *acc, size_t len,
                double *restrict out)
{
    const double *evidence  = accum_const_data(acc, ACCUM_SPANNED);
    const double *mutations = accum_const_data(acc, ACCUM_MUTATIONS);
    const double *squared   = accum_const_data(acc, ACCUM_MUTATIONS_SQUARED);

    for (size_t i = 0; i < len; i++) {
        out[i] = error_at(cfg, mutations[i], squared[i], evidence[i], i, len);
    }
}
