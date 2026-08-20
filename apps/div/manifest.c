/* manifest.c -- the datasets cmuts div writes, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* The datasets a run leaves behind, which are those every output holds: what it was
 * given, divided by the control it was given. */
static const char SUMMED[] = "Summed over the inputs.";

static const fmt_written WRITTEN[] = {
    { .id = FMT_REACTIVITY,
      .note = "Sample over control.",
      .depends = FMT_DEPENDS(FMT_REACTIVITY),
      .required = true },
    { .id = FMT_ERROR,
      .note = "Propagated from the inputs.",
      .depends = FMT_DEPENDS(FMT_REACTIVITY, FMT_ERROR) },
    { .id = FMT_COVERAGE,  .note = SUMMED, .depends = FMT_DEPENDS(FMT_COVERAGE) },
    { .id = FMT_SEQUENCE,  .depends = FMT_DEPENDS(FMT_SEQUENCE) },
    { .id = FMT_LENGTHS,   .note = SUMMED, .depends = FMT_DEPENDS(FMT_LENGTHS) },
    { .id = FMT_READS,     .note = SUMMED, .depends = FMT_DEPENDS(FMT_READS) },
    { .id = FMT_REJECTED,  .note = SUMMED, .depends = FMT_DEPENDS(FMT_REJECTED) },
    { .id = FMT_UNMAPPED,  .note = SUMMED, .depends = FMT_DEPENDS(FMT_UNMAPPED) },
};

const fmt_manifest CMUTS_DIV_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };
