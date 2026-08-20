/* manifest.c -- the datasets cmuts hmm writes, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* The datasets a run leaves behind. Every output holds the first of them; the pairwise
 * squares are written only when they are asked for. */
static const out_written WRITTEN[] = {
    { .id = OUT_REACTIVITY },
    { .id = OUT_ERROR },
    { .id = OUT_COVERAGE },
    { .id = OUT_SEQUENCE },
    { .id = OUT_LENGTHS },
    { .id = OUT_READS },
    { .id = OUT_REJECTED },
    { .id = OUT_UNMAPPED },
    { .id = OUT_PAIRWISE_CORRELATION, .condition = "--pairwise correlation" },
    { .id = OUT_PAIRWISE_CONDITIONAL, .condition = "--pairwise conditional" },
    { .id = OUT_PAIRWISE_COVERAGE,    .condition = "--pairwise" },
};


const out_manifest CMUTS_HMM_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };
