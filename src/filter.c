/* filter.c -- alignment acceptance tests.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "filter.h"

#include "align.h"

/* The mapping quality a read must reach by default. A MAPQ is the aligner's phred-scaled
 * chance that a placement is wrong, so this admits a read placed with better than one
 * chance in a hundred of being placed elsewhere. */
#define DEFAULT_MIN_MAPQ 20

filter_config filter_defaults(void)
{
    return (filter_config){
        .min_mapq   = DEFAULT_MIN_MAPQ,
        .strand     = FILTER_STRAND_FORWARD | FILTER_STRAND_REVERSE,
        .min_length = FILTER_LENGTH_UNBOUNDED,
        .max_length = FILTER_LENGTH_UNBOUNDED,
    };
}

bool filter_satisfiable(const filter_config *filter)
{
    return filter->max_length == FILTER_LENGTH_UNBOUNDED
           || filter->max_length >= filter->min_length;
}

/* MAPQ 255 means "unavailable" in the SAM specification, so a read carrying it is discarded
 * at every threshold: an aligner that reported no confidence in a placement is not reporting
 * confidence of 255. samtools compares it numerically, admitting it everywhere, so the two
 * disagree on such a read. */
static bool mapping_quality_accepted(const filter_config *filter, const cm_bam_record *read)
{
    return read->mapq != FILTER_MAPQ_UNAVAILABLE && read->mapq >= filter->min_mapq;
}

static bool strand_accepted(const filter_config *filter, const cm_bam_record *read)
{
    bool reverse = (read->flag & BAM_FREVERSE) != 0;
    int  own     = reverse ? FILTER_STRAND_REVERSE : FILTER_STRAND_FORWARD;

    return (filter->strand & own) != 0;
}

/* Returns whether the read falls inside the length bounds. The length is that of the
 * stored sequence, so a hard-clipped read counts only the bases the aligner kept. */
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

/* Returns whether the record stores a sequence. SAM spells its absence as a SEQ and
 * QUAL of "*", leaving no bases to compare against the reference. */
static bool sequence_present(const cm_bam_record *read)
{
    return read->seq != NULL;
}

/* Returns whether the record places any of its read. It may place none three ways: no
 * CIGAR at all, a CIGAR made wholly of clips, or one consuming reference but no read.
 * There is then no position for the read to contribute to. */
static bool placement_present(const cm_bam_record *read)
{
    aln_span placed = aln_placed_span(read);

    return placed.end > placed.begin;
}

/* Returns whether the record is a primary alignment. A secondary one places a read
 * already counted at its primary, so accepting it would count one molecule twice.
 * Supplementary alignments carry distinct pieces of a split read and are not excluded. */
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
