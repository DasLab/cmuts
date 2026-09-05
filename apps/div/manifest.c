/* manifest.c -- the datasets cmuts div writes, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* The datasets a run leaves behind, which are those every output holds: what it was
 * given, divided by the control it was given. */
static const char SUMMED[] = "Summed over the inputs.";

static const fmt_written WRITTEN[] = {
    { .id = FMT_MISMATCH_RATE,
      .how = "Sample over control.",
      .depends = FMT_DEPENDS(FMT_MISMATCH_RATE),
      .required = true },
    { .id = FMT_MISMATCH_ERROR,
      .how = "Propagated from the inputs.",
      .depends = FMT_DEPENDS(FMT_MISMATCH_RATE, FMT_MISMATCH_ERROR) },
    { .id = FMT_INSERTION_RATE,
      .how = "Sample over control.",
      .depends = FMT_DEPENDS(FMT_INSERTION_RATE),
      .required = true },
    { .id = FMT_INSERTION_ERROR,
      .how = "Propagated from the inputs.",
      .depends = FMT_DEPENDS(FMT_INSERTION_RATE, FMT_INSERTION_ERROR) },
    { .id = FMT_DELETION_RATE,
      .how = "Sample over control.",
      .depends = FMT_DEPENDS(FMT_DELETION_RATE),
      .required = true },
    { .id = FMT_DELETION_ERROR,
      .how = "Propagated from the inputs.",
      .depends = FMT_DEPENDS(FMT_DELETION_RATE, FMT_DELETION_ERROR) },
    { .id = FMT_COVERAGE,  .how = SUMMED, .depends = FMT_DEPENDS(FMT_COVERAGE) },
    { .id = FMT_SEQUENCE,
      .how = "Copied from the inputs.",
      .depends = FMT_DEPENDS(FMT_SEQUENCE) },
    { .id = FMT_LENGTHS,   .how = SUMMED, .depends = FMT_DEPENDS(FMT_LENGTHS) },
    { .id = FMT_READS,     .how = SUMMED, .depends = FMT_DEPENDS(FMT_READS) },
    { .id = FMT_REJECTED,  .how = SUMMED, .depends = FMT_DEPENDS(FMT_REJECTED) },
    { .id = FMT_UNMAPPED,  .how = SUMMED, .depends = FMT_DEPENDS(FMT_UNMAPPED) },
};

const fmt_manifest CMUTS_DIV_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };
