/* manifest.c -- the datasets cmuts norm writes, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* The datasets a run leaves behind: what it was given, divided by one norm, and that
 * norm alongside it. */
static const char DIVIDED[] = "Divided by the norm.";
static const char COPIED[]  = "Copied from the input.";

static const fmt_written WRITTEN[] = {
    { .id = FMT_MISMATCH_RATE,
      .how = DIVIDED,
      .depends = FMT_DEPENDS(FMT_MISMATCH_RATE),
      .required = true },
    { .id = FMT_MISMATCH_ERROR,
      .how = DIVIDED,
      .depends = FMT_DEPENDS(FMT_MISMATCH_ERROR) },
    { .id = FMT_NORM,
      .how = "Estimated per the specified scheme.",
      .depends = FMT_DEPENDS(FMT_MISMATCH_RATE, FMT_COVERAGE),
      .required = true },
    { .id = FMT_COVERAGE,  .how = COPIED, .depends = FMT_DEPENDS(FMT_COVERAGE) },
    { .id = FMT_SEQUENCE,  .how = COPIED, .depends = FMT_DEPENDS(FMT_SEQUENCE) },
    { .id = FMT_LENGTHS,   .how = COPIED, .depends = FMT_DEPENDS(FMT_LENGTHS) },
    { .id = FMT_READS,     .how = COPIED, .depends = FMT_DEPENDS(FMT_READS) },
    { .id = FMT_REJECTED,  .how = COPIED, .depends = FMT_DEPENDS(FMT_REJECTED) },
    { .id = FMT_UNMAPPED,  .how = COPIED, .depends = FMT_DEPENDS(FMT_UNMAPPED) },
};

const fmt_manifest CMUTS_NORM_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };
