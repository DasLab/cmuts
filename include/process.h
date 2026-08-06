/* process.h -- per-read contribution to a reference's accumulated values.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "accum.h"
#include "bam.h"
#include "fasta.h"
#include "phred.h"

/* What the processing step is told besides the read itself: how to read a
 * quality, and in time what to count as a modification and what to weigh it by.
 *
 * Built before any worker starts and never written afterwards, which is what
 * lets every worker be handed the same one. */
typedef struct {
    phred quality;
} process_config;

void process_config_build(process_config *config);

/* Adds one read's contribution to target, which is never cleared here.
 *
 * Accumulating into a caller-provided target rather than returning a fresh
 * struct is what lets a worker hold a private accumulator across a run of
 * reads: contributions land there with no allocation and no locking, and reach
 * the shared per-reference accumulator in a single merge. */
void process(const cm_bam_record *read, const cm_fasta_record *ref,
             const process_config *config, accum *target);
