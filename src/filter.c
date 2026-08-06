/* filter.c -- alignment acceptance tests.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "filter.h"

filter_config filter_defaults(void)
{
    return (filter_config){
        .min_mapq   = 0,
        .strand     = FILTER_STRAND_BOTH,
        .min_length = FILTER_LENGTH_UNBOUNDED,
        .max_length = FILTER_LENGTH_UNBOUNDED,
    };
}

/* MAPQ 255 means "unavailable" in the SAM specification rather than "perfect",
 * but it is compared numerically here, as samtools does. An aligner that emits
 * 255 throughout would otherwise have the whole of its output discarded by any
 * threshold at all. */
static bool mapping_quality_accepted(const filter_config *filter, const cm_bam_record *read)
{
    return read->mapq >= filter->min_mapq;
}

static bool strand_accepted(const filter_config *filter, const cm_bam_record *read)
{
    bool reverse = (read->flag & BAM_FREVERSE) != 0;

    switch (filter->strand) {
        case FILTER_STRAND_FORWARD: return !reverse;
        case FILTER_STRAND_REVERSE: return reverse;
        default:                    return true;
    }
}

/* Length is that of the stored sequence, so a hard-clipped read counts only the
 * bases the aligner kept, not those it trimmed away. */
static bool length_accepted(const filter_config *filter, const cm_bam_record *read)
{
    if (filter->min_length != FILTER_LENGTH_UNBOUNDED && read->l_qseq < filter->min_length)
        return false;

    if (filter->max_length != FILTER_LENGTH_UNBOUNDED && read->l_qseq > filter->max_length)
        return false;

    return true;
}

bool filter_accepts(const filter_config *filter, const cm_bam_record *read)
{
    return mapping_quality_accepted(filter, read) &&
           strand_accepted(filter, read) &&
           length_accepted(filter, read);
}
