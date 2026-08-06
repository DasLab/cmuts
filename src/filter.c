/* filter.c -- alignment acceptance tests.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "filter.h"

filter_config filter_defaults(void)
{
    return (filter_config){ .min_mapq = 0 };
}

/* MAPQ 255 means "unavailable" in the SAM specification rather than "perfect",
 * but it is compared numerically here, as samtools does. An aligner that emits
 * 255 throughout would otherwise have the whole of its output discarded by any
 * threshold at all. */
bool filter_accepts(const filter_config *filter, const cm_bam_record *read)
{
    return read->mapq >= filter->min_mapq;
}
