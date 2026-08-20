/* manifest.c -- the datasets cmuts score reads, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

static const out_written READS[] = {
    { .id = OUT_COVERAGE,   .origin = OUT_REQUIRED },
    { .id = OUT_REACTIVITY, .origin = OUT_REQUIRED },
};

const out_manifest CMUTS_SCORE_READS = { READS, sizeof READS / sizeof *READS };
