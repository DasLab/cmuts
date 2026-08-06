/* filter.h -- which alignments reach the processing step.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>

#include "bam.h"

/* Criteria an alignment must meet to be processed.
 *
 * Unmapped reads are excluded before any of this and counted separately, since
 * they belong to no reference and so have nowhere to be accumulated. */
typedef struct {
    int min_mapq;  /* 0 to 255; alignments scoring below it are discarded */
} filter_config;

filter_config filter_defaults(void);

bool filter_accepts(const filter_config *filter, const cm_bam_record *read);
