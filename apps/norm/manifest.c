/* manifest.c -- the datasets cmuts norm writes, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* The datasets a run leaves behind: what it was given, divided by one norm, and that
 * norm alongside it. */
static const fmt_written WRITTEN[] = {
    { .id = FMT_REACTIVITY,
      .note = "Divided by the norm.",
      .depends = FMT_DEPENDS(FMT_REACTIVITY),
      .required = true },
    { .id = FMT_ERROR,
      .note = "Divided by the norm.",
      .depends = FMT_DEPENDS(FMT_ERROR) },
    { .id = FMT_NORM,
      .depends = FMT_DEPENDS(FMT_REACTIVITY, FMT_COVERAGE),
      .required = true },
    { .id = FMT_COVERAGE,  .depends = FMT_DEPENDS(FMT_COVERAGE) },
    { .id = FMT_SEQUENCE,  .depends = FMT_DEPENDS(FMT_SEQUENCE) },
    { .id = FMT_LENGTHS,   .depends = FMT_DEPENDS(FMT_LENGTHS) },
    { .id = FMT_READS,     .depends = FMT_DEPENDS(FMT_READS) },
    { .id = FMT_REJECTED,  .depends = FMT_DEPENDS(FMT_REJECTED) },
    { .id = FMT_UNMAPPED,  .depends = FMT_DEPENDS(FMT_UNMAPPED) },
};

const fmt_manifest CMUTS_NORM_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };
