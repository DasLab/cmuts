/* manifest.c -- the datasets cmuts hmm writes, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* The datasets a run leaves behind. Every output holds the first of them; the pairwise
 * squares are written only when they are asked for. */
static const fmt_written WRITTEN[] = {
    { .id = FMT_REACTIVITY },
    { .id = FMT_ERROR },
    { .id = FMT_COVERAGE },
    { .id = FMT_SEQUENCE },
    { .id = FMT_LENGTHS },
    { .id = FMT_READS },
    { .id = FMT_REJECTED },
    { .id = FMT_UNMAPPED },
    { .id = FMT_PAIRWISE_CORRELATION, .condition = "--pairwise correlation" },
    { .id = FMT_PAIRWISE_CONDITIONAL, .condition = "--pairwise conditional" },
    { .id = FMT_PAIRWISE_COVERAGE,    .condition = "--pairwise" },
};


const fmt_manifest CMUTS_HMM_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };
