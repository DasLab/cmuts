/* manifest.c -- the datasets cmuts score reads, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

const fmt_reads CMUTS_SCORE_READS = {
    .fields = {
        { FMT_COVERAGE,   .required = true },
        { FMT_REACTIVITY, .required = true },
    },
    .n_fields = 2,
};
