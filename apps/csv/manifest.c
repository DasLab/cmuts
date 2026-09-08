/* manifest.c -- the datasets cmuts csv reads, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* Every per-base field of an output, in the order format.h lists them, which is the order
 * of the columns. A channel that the run did not write is omitted from the table, and does
 * not cause the input to be refused.
 *
 * The coverage is required, because every output holds it. The sequence is required
 * because it supplies the base at each position and the length of each reference. */
const fmt_request CMUTS_CSV_READS[] = {
    { .id = FMT_MISMATCH_RATE   },
    { .id = FMT_MISMATCH_ERROR  },
    { .id = FMT_INSERTION_RATE  },
    { .id = FMT_INSERTION_ERROR },
    { .id = FMT_DELETION_RATE   },
    { .id = FMT_DELETION_ERROR  },
    { .id = FMT_COVERAGE, .required = true },
    { .id = FMT_SEQUENCE, .required = true },
};

const size_t CMUTS_CSV_N_READS = sizeof CMUTS_CSV_READS / sizeof *CMUTS_CSV_READS;
