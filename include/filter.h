/* filter.h -- which alignments reach the processing step.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>

#include "bam.h"

/* Which strands' alignments to keep, as the set of them. The test is on the alignment's
 * own reverse bit, which describes the read and not the fragment. */
typedef enum {
    FILTER_STRAND_FORWARD = 1 << 0,
    FILTER_STRAND_REVERSE = 1 << 1,
} filter_strand;

/* Each bit selects an alignment type to keep. A supplementary record holds part
 * of a split read, and the primary record holds the rest. */
typedef enum {
    FILTER_ALIGNMENT_PRIMARY       = 1 << 0,
    FILTER_ALIGNMENT_SUPPLEMENTARY = 1 << 1,
} filter_alignment;

/* A length bound of this is not applied at all. */
#define FILTER_LENGTH_UNBOUNDED 0

/* The mapping quality of a placement the aligner reported no confidence in. */
#define FILTER_MAPQ_UNAVAILABLE 255

/* Criteria an alignment must meet to be processed.
 *
 * Secondary alignments, records storing no sequence, records carrying no CIGAR
 * and placements of unavailable mapping quality are always excluded. */
typedef struct {
    int min_mapq;       /* 0 to 254; alignments scoring below it are discarded */
    int strand;         /* filter_strand bits; a read on neither strand is discarded */
    int alignment_type; /* Records of a type not in these bits are discarded. */
    int min_length;     /* FILTER_LENGTH_UNBOUNDED for no lower bound */
    int max_length;     /* FILTER_LENGTH_UNBOUNDED for no upper bound */
} filter_config;

filter_config filter_defaults(void);

/* Returns whether any read can meet the criteria: an upper length bound below the lower
 * one is met by no read of any length. */
bool filter_satisfiable(const filter_config *filter);

bool filter_accepts(const filter_config *filter, const cm_bam_record *read);

/* Returns whether the record places a further piece of a read that another record already
 * places. A total over reads counts the read once, so it passes over these. */
bool filter_is_supplementary(const cm_bam_record *read);
