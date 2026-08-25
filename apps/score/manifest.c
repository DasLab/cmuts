/* manifest.c -- the datasets cmuts score reads, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

const fmt_request CMUTS_SCORE_READS[] = {
    { FMT_COVERAGE,   .required = true },
    { FMT_MISMATCH_RATE, .required = true },
};

const size_t CMUTS_SCORE_N_READS = sizeof CMUTS_SCORE_READS / sizeof *CMUTS_SCORE_READS;
