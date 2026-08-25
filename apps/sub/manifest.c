/* manifest.c -- the datasets cmuts sub writes, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* The datasets a run leaves behind, which are those every output holds: what it was
 * given, less the background it was given. */
static const char SUMMED[] = "Summed over the inputs.";

static const fmt_written WRITTEN[] = {
    { .id = FMT_MISMATCH_RATE,
      .how = "Treated less untreated.",
      .depends = FMT_DEPENDS(FMT_MISMATCH_RATE),
      .required = true },
    { .id = FMT_MISMATCH_ERROR,
      .how = "Propagated from the inputs.",
      .depends = FMT_DEPENDS(FMT_MISMATCH_ERROR) },
    { .id = FMT_COVERAGE,  .how = SUMMED, .depends = FMT_DEPENDS(FMT_COVERAGE) },
    { .id = FMT_SEQUENCE,
      .how = "Copied from the inputs.",
      .depends = FMT_DEPENDS(FMT_SEQUENCE) },
    { .id = FMT_LENGTHS,   .how = SUMMED, .depends = FMT_DEPENDS(FMT_LENGTHS) },
    { .id = FMT_READS,     .how = SUMMED, .depends = FMT_DEPENDS(FMT_READS) },
    { .id = FMT_REJECTED,  .how = SUMMED, .depends = FMT_DEPENDS(FMT_REJECTED) },
    { .id = FMT_UNMAPPED,  .how = SUMMED, .depends = FMT_DEPENDS(FMT_UNMAPPED) },
};

const fmt_manifest CMUTS_SUB_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };
