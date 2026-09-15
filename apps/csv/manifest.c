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
#define CHANNEL(channel, rate, error) \
    { .id = (rate) },                 \
    { .id = (error) },

const fmt_request CMUTS_CSV_READS[] = {
    FMT_CHANNELS(CHANNEL)
    { .id = FMT_COVERAGE, .required = true },
    { .id = FMT_SEQUENCE, .required = true },
};

#undef CHANNEL

const size_t CMUTS_CSV_N_READS = sizeof CMUTS_CSV_READS / sizeof *CMUTS_CSV_READS;
