/* manifest.c -- the datasets cmuts div writes, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* The datasets a run leaves behind, which are those every output holds: what it was
 * given, divided by the control it was given. */
static const char SUMMED[] = "Summed over the inputs.";

static const out_written WRITTEN[] = {
    { .id = OUT_REACTIVITY,
      .note = "Sample over control.",
      .origin = OUT_REQUIRED },
    { .id = OUT_ERROR,
      .note = "Propagated from the inputs." },
    { .id = OUT_COVERAGE,  .note = SUMMED },
    { .id = OUT_SEQUENCE },
    { .id = OUT_LENGTHS,   .note = SUMMED },
    { .id = OUT_READS,     .note = SUMMED },
    { .id = OUT_REJECTED,  .note = SUMMED },
    { .id = OUT_UNMAPPED,  .note = SUMMED },
};


const out_manifest CMUTS_DIV_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };
