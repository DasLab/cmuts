/* filter.c -- alignment acceptance tests.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "filter.h"

#include "align.h"

filter_config filter_defaults(void)
{
    return (filter_config){
        .min_mapq   = 20,
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

/* A record may store no sequence at all, which SAM spells as a QUAL and SEQ of
 * "*". There is nothing to compare against the reference, so such a read is
 * turned away here rather than reaching a step that would have to invent an
 * answer for it. */
static bool sequence_present(const cm_bam_record *read)
{
    return read->seq != NULL;
}

/* A record may place none of its read: no CIGAR at all, which SAM spells the
 * same way as an absent sequence; one made wholly of clips, naming bases the
 * aligner declined to put anywhere; or one naming no base of the read at all,
 * having nothing but reference in it. There is then no position for the read to
 * say anything about, so it is turned away here rather than counted as one that
 * went on to contribute nothing to any of them. */
static bool placement_present(const cm_bam_record *read)
{
    aln_span placed = aln_placed_span(read);

    return placed.end > placed.begin;
}

/* A secondary alignment places a read already counted at its primary, so
 * accepting one would count a single molecule twice. Supplementary alignments
 * carry distinct pieces of a split read, and are not refused here. */
static bool is_primary(const cm_bam_record *read)
{
    return (read->flag & BAM_FSECONDARY) == 0;
}

bool filter_accepts(const filter_config *filter, const cm_bam_record *read)
{
    return sequence_present(read) &&
           placement_present(read) &&
           is_primary(read) &&
           mapping_quality_accepted(filter, read) &&
           strand_accepted(filter, read) &&
           length_accepted(filter, read);
}
