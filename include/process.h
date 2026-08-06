/* process.h -- per-read contribution to a reference's accumulated values.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "accum.h"
#include "bam.h"
#include "fasta.h"

/* Adds one read's contribution to target, which is never cleared here.
 *
 * Accumulating into a caller-provided target rather than returning a fresh
 * struct is what lets a worker hold a private accumulator across a run of
 * reads: contributions land there with no allocation and no locking, and reach
 * the shared per-reference accumulator in a single merge. */
void process(const cm_bam_record *read, const cm_fasta_record *ref, accum *target);
