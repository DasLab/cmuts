/* manifest.c -- the datasets cmuts hmm writes, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* The datasets a run leaves behind. Every output holds the first of them; the pairwise
 * squares are written only when they are asked for. */
static const char HMM[] = "Estimated by the HMM.";

static const char BINOMIAL[] =
    "The binomial standard error of the rate at the position's depth.";

/* Every channel is estimated and reported the same way. */
#define CHANNEL(channel, rate, error) \
    { .id = (rate),  .how = HMM },      \
    { .id = (error), .how = BINOMIAL },

static const fmt_written WRITTEN[] = {
    FMT_CHANNELS(CHANNEL)
    { .id = FMT_COVERAGE,  .how = HMM },
    { .id = FMT_SEQUENCE,  .how = "Tokenized from the FASTA." },
    { .id = FMT_LENGTHS,   .how = "The read length reported in the alignment." },
    { .id = FMT_READS,     .how = "The number of reads the HMM successfully processed." },
    { .id = FMT_REJECTED,
      .how = "The number of reads rejected by a filter or by the HMM." },
    { .id = FMT_UNMAPPED,
      .how = "The number of unmapped reads in the alignment." },
    { .id = FMT_PAIRWISE_CORRELATION, .how = HMM, .condition = "--pairwise correlation" },
    { .id = FMT_PAIRWISE_CONDITIONAL, .how = HMM, .condition = "--pairwise conditional" },
    { .id = FMT_PAIRWISE_COVERAGE,
      .how = HMM,
      .condition = "--pairwise" },
};


#undef CHANNEL

const fmt_manifest CMUTS_HMM_WRITES = {
    WRITTEN, sizeof WRITTEN / sizeof *WRITTEN, NULL,
};
