/* manifest.c -- the datasets cmuts norm writes, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "manifest.h"

/* The datasets a run leaves behind: what it was given, divided by one norm, and that
 * norm alongside it. */
static const out_written WRITTEN[] = {
    { .id = OUT_REACTIVITY,
      .note = "Divided by the norm.",
      .origin = OUT_REQUIRED },
    { .id = OUT_ERROR,
      .note = "Divided by the norm." },
    { .id = OUT_NORM,
      .origin = OUT_MADE },
    { .id = OUT_COVERAGE,
      .origin = OUT_REQUIRED },
    { .id = OUT_SEQUENCE },
    { .id = OUT_LENGTHS },
    { .id = OUT_READS },
    { .id = OUT_REJECTED },
    { .id = OUT_UNMAPPED },
};


const out_manifest CMUTS_NORM_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };
