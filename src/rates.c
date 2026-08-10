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

/* Some evidence is wanted whatever depth was asked for, so a depth of nothing
 * means whatever there is rather than none at all. */
static bool known_at(double wanted, double evidence)
{
    return evidence > 0.0 && evidence >= wanted;
}

/* The rate cannot exceed one: every weight is a share of an event, and an
 * insertion spans what it lays. It is held there all the same, so that what
 * comes of it is a proportion however the last bit of a double rounded. */
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

/* The standard error of a proportion over the evidence standing as its count. */
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
