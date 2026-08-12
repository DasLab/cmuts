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

/* MAPQ 255 means "unavailable" in the SAM specification rather than "perfect", but is
 * compared numerically here as samtools does. An aligner emitting 255 throughout would
 * otherwise have all of its output discarded by any threshold at all. */
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

/* Length is that of the stored sequence, so a hard-clipped read counts only the bases
 * the aligner kept. */
static bool length_accepted(const filter_config *filter, const cm_bam_record *read)
{
    if (filter->min_length != FILTER_LENGTH_UNBOUNDED && read->l_qseq < filter->min_length) {
        return false;
    }

    if (filter->max_length != FILTER_LENGTH_UNBOUNDED && read->l_qseq > filter->max_length) {
        return false;
    }

    return true;
}

/* Whether the record stores a sequence. SAM spells its absence as a SEQ and QUAL of
 * "*", leaving nothing to compare against the reference. */
static bool sequence_present(const cm_bam_record *read)
{
    return read->seq != NULL;
}

/* Whether the record places any of its read. It may place none three ways: no CIGAR at
 * all, a CIGAR made wholly of clips, or one consuming reference but no read. There is
 * then no position for the read to contribute to. */
static bool placement_present(const cm_bam_record *read)
{
    aln_span placed = aln_placed_span(read);

    return placed.end > placed.begin;
}

/* Whether the record is a primary alignment. A secondary one places a read already
 * counted at its primary, so accepting it would count one molecule twice. Supplementary
 * alignments carry distinct pieces of a split read and are not excluded. */
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
