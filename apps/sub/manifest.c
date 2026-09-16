/* manifest.c -- the datasets cmuts sub writes, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* The datasets a run leaves behind, which are those every output holds: what it was
 * given, less the background it was given. */
static const char SUMMED[]     = "Summed over the inputs.";
static const char SUBTRACTED[] = "Treated less untreated.";
static const char PROPAGATED[] = "Propagated from the inputs.";

/* Every channel is subtracted the same way. Each field is written only where every
 * input holds it. */
#define CHANNEL(channel, rate, error) \
    { .id = (rate),                   \
      .how = SUBTRACTED,              \
      .depends = FMT_DEPENDS(rate) }, \
    { .id = (error),                  \
      .how = PROPAGATED,              \
      .depends = FMT_DEPENDS(error) },

/* At least one rate is required. */
#define ANY_RATE(channel, rate, error) rate,

static const fmt_field_id ANY_OF[] = { FMT_CHANNELS(ANY_RATE) FMT_N_FIELDS };

#undef ANY_RATE

static const fmt_written WRITTEN[] = {
    FMT_CHANNELS(CHANNEL)
    { .id = FMT_COVERAGE,  .how = SUMMED, .depends = FMT_DEPENDS(FMT_COVERAGE) },
    { .id = FMT_SEQUENCE,
      .how = "Copied from the inputs.",
      .depends = FMT_DEPENDS(FMT_SEQUENCE) },
    { .id = FMT_LENGTHS,   .how = SUMMED, .depends = FMT_DEPENDS(FMT_LENGTHS) },
    { .id = FMT_PRIMARY_COUNTED,
      .how = SUMMED, .depends = FMT_DEPENDS(FMT_PRIMARY_COUNTED) },
    { .id = FMT_PRIMARY_REJECTED,
      .how = SUMMED, .depends = FMT_DEPENDS(FMT_PRIMARY_REJECTED) },
    { .id = FMT_SUPPLEMENTARY_COUNTED,
      .how = SUMMED, .depends = FMT_DEPENDS(FMT_SUPPLEMENTARY_COUNTED) },
    { .id = FMT_SUPPLEMENTARY_REJECTED,
      .how = SUMMED, .depends = FMT_DEPENDS(FMT_SUPPLEMENTARY_REJECTED) },
    { .id = FMT_UNMAPPED,  .how = SUMMED, .depends = FMT_DEPENDS(FMT_UNMAPPED) },
};

#undef CHANNEL

const fmt_manifest CMUTS_SUB_WRITES = {
    WRITTEN, sizeof WRITTEN / sizeof *WRITTEN, ANY_OF,
};
