/* filter.h -- which alignments reach the processing step.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>

#include "bam.h"

/* Which strand's alignments to keep.
 *
 * The test is on the alignment's own reverse bit, which for single-end reads
 * is the strand the read came from. It is not a statement about the fragment:
 * for paired data, what strand a fragment belongs to depends on the library
 * protocol and on which mate is in hand. */
typedef enum {
    FILTER_STRAND_BOTH,
    FILTER_STRAND_FORWARD,
    FILTER_STRAND_REVERSE,
} filter_strand;

/* Criteria an alignment must meet to be processed.
 *
 * Unmapped reads are excluded before any of this and counted separately, since
 * they belong to no reference and so have nowhere to be accumulated. Secondary
 * alignments and records storing no sequence are refused whatever is set here.
 *
 * The fields are int rather than narrower or enumerated types because the
 * command line writes them directly, through a pointer of the declared type. */
/* A length bound of this is not applied at all. It doubles as the identity for
 * either bound: no read is shorter than zero, and an absent upper bound is not
 * a length. */
#define FILTER_LENGTH_UNBOUNDED 0

typedef struct {
    int min_mapq;    /* 0 to 255; alignments scoring below it are discarded */
    int strand;      /* a filter_strand */
    int min_length;  /* FILTER_LENGTH_UNBOUNDED for no lower bound */
    int max_length;  /* FILTER_LENGTH_UNBOUNDED for no upper bound */
} filter_config;

filter_config filter_defaults(void);

bool filter_accepts(const filter_config *filter, const cm_bam_record *read);
