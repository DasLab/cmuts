/* manifest.c -- the datasets cmuts norm writes, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* The datasets a run leaves behind: what it was given, divided by one norm, and that
 * norm alongside it. */
static const char DIVIDED[] = "Divided by the norm.";
static const char COPIED[]  = "Copied from the input.";

/* Every channel takes the norm. Each field is written only where the input holds it. */
#define CHANNEL(channel, rate, error) \
    { .id = (rate),                   \
      .how = DIVIDED,                 \
      .depends = FMT_DEPENDS(rate) }, \
    { .id = (error),                  \
      .how = DIVIDED,                 \
      .depends = FMT_DEPENDS(error) },

/* At least one of these rates is required. The aggregate is taken over whichever of them
 * the input holds. */
static const fmt_field_id ANY_OF[] = {
    FMT_MISMATCH_RATE, FMT_INSERTION_RATE, FMT_DELETION_RATE, FMT_N_FIELDS,
};

static const fmt_written WRITTEN[] = {
    FMT_CHANNELS(CHANNEL)
    { .id = FMT_NORM,
      .how = "Estimated per the specified scheme, over the aggregate of whichever mismatch, insertion and deletion rates the input holds.",
      .depends = FMT_DEPENDS(FMT_COVERAGE),
      .required = true },
    { .id = FMT_COVERAGE,  .how = COPIED, .depends = FMT_DEPENDS(FMT_COVERAGE) },
    { .id = FMT_SEQUENCE,  .how = COPIED, .depends = FMT_DEPENDS(FMT_SEQUENCE) },
    { .id = FMT_LENGTHS,   .how = COPIED, .depends = FMT_DEPENDS(FMT_LENGTHS) },
    { .id = FMT_READS,     .how = COPIED, .depends = FMT_DEPENDS(FMT_READS) },
    { .id = FMT_REJECTED,  .how = COPIED, .depends = FMT_DEPENDS(FMT_REJECTED) },
    { .id = FMT_UNMAPPED,  .how = COPIED, .depends = FMT_DEPENDS(FMT_UNMAPPED) },
};

#undef CHANNEL

const fmt_manifest CMUTS_NORM_WRITES = {
    WRITTEN, sizeof WRITTEN / sizeof *WRITTEN, ANY_OF,
};
