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

/* Whether a position carries enough evidence to report on. Some is required whatever
 * depth was asked for, so a depth of zero means whatever there is and not none. */
static bool known_at(double wanted, double evidence)
{
    return evidence > 0.0 && evidence >= wanted;
}

/* The mutations at a position over the evidence for them, held to one. Every weight is
 * a share of an event and an insertion spans what it contributes, so the ratio cannot
 * exceed one except by rounding. */
static double rate_of(double mutations, double evidence)
{
    double rate = evidence > 0.0 ? mutations / evidence : 0.0;

    return rate > 1.0 ? 1.0 : rate;
}

void rate_reactivity(const rate_config *cfg, const accum *acc, size_t len,
                     double *restrict out)
{
    const double *evidence  = accum_const_data(acc, ACCUM_SPANNED);
    const double *mutations = accum_const_data(acc, ACCUM_MUTATIONS);
    double        wanted    = cfg->min_depth;

    for (size_t i = 0; i < len; i++)
        out[i] = known_at(wanted, evidence[i])
               ? rate_of(mutations[i], evidence[i])
               : (double)NAN;
}

/* The standard error of the rate, taking the evidence as the count it is a proportion
 * of. */
void rate_error(const rate_config *cfg, const accum *acc, size_t len,
                double *restrict out)
{
    const double *evidence  = accum_const_data(acc, ACCUM_SPANNED);
    const double *mutations = accum_const_data(acc, ACCUM_MUTATIONS);
    double        wanted    = cfg->min_depth;

    for (size_t i = 0; i < len; i++) {
        double rate = rate_of(mutations[i], evidence[i]);

        out[i] = known_at(wanted, evidence[i])
               ? sqrt(rate * (1.0 - rate) / evidence[i])
               : (double)NAN;
    }
}
