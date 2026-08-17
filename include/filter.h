/* filter.h -- which alignments reach the processing step.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>

#include "bam.h"

/* Which strands' alignments to keep, as the set of them: keeping both is asking for each.
 *
 * The test is on the alignment's own reverse bit, which for single-end reads is the strand
 * the read came from. It describes the read and not the fragment: for paired data, what
 * fragment belongs to depends on the library protocol and on which mate is being read. */
typedef enum {
    FILTER_STRAND_FORWARD = 1 << 0,
    FILTER_STRAND_REVERSE = 1 << 1,
} filter_strand;

/* A length bound of this is not applied at all. It also serves as the identity for either
 * bound: no read is shorter than zero, and an absent upper bound is not a length. */
#define FILTER_LENGTH_UNBOUNDED 0

/* The mapping quality of a placement the aligner reported no confidence in. */
#define FILTER_MAPQ_UNAVAILABLE 255

/* Criteria an alignment must meet to be processed.
 *
 * Unmapped reads are excluded before any of this and counted separately, belonging to no
 * reference and so having nowhere to be accumulated. Secondary alignments, records storing
 * no sequence, records carrying no CIGAR and placements of unavailable mapping quality are
 * excluded whatever is set here.
 *
 * The fields are int because the command line writes them directly, through a pointer of
 * the declared type. */
typedef struct {
    int min_mapq;    /* 0 to 254; alignments scoring below it are discarded */
    int strand;      /* filter_strand bits; a read on neither strand is discarded */
    int min_length;  /* FILTER_LENGTH_UNBOUNDED for no lower bound */
    int max_length;  /* FILTER_LENGTH_UNBOUNDED for no upper bound */
} filter_config;

filter_config filter_defaults(void);

/* Returns whether any read can meet the criteria. An upper length bound below the lower one
 * is met by no read of any length, unlike a bound that merely no read in a file happens to
 * meet. */
bool filter_satisfiable(const filter_config *filter);

bool filter_accepts(const filter_config *filter, const cm_bam_record *read);
