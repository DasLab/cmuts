/* tally.h -- one read's contribution to a reference's accumulated values.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "accum.h"
#include "bam.h"
#include "fasta.h"
#include "pairs.h"
#include "phmm.h"
#include "phred.h"

/* The settings a caller supplies, in types a command line can write into directly. */
typedef struct {
    int          band;     /* reference positions the marginal may look either side of the
                              CIGAR; 0 pins it to the path as written */
    phmm_weights weights;  /* what an event of each kind counts for */
    phmm_params  params;   /* the rates the model runs on */
} tally_config;

tally_config tally_defaults(void);

/* What the tally works from besides the read itself: the quality table and the model.
 *
 * Derived from the config, built before any worker starts and never written afterwards,
 * so every worker may be given the same one. */
typedef struct {
    phred quality;
    phmm  model;
    int   band;   /* the half-width every row of the marginal is given */
} tally_tables;

void tally_tables_build(tally_tables *tables, const tally_config *config);

/* Working buffers one worker reuses across every read it processes.
 *
 * Private to a worker, since the marginal writes over the whole of it for each read, and
 * reused across reads. Every read is counted through one, so a caller must create it
 * before any read arrives. */
typedef struct tally_scratch tally_scratch;

tally_scratch *tally_scratch_create(void);
void           tally_scratch_destroy(tally_scratch *scratch);

/* Adds one read's contribution to target, which is never cleared here.
 *
 * A caller-provided target lets a worker hold a private accumulator across a run of
 * reads: contributions land there with no allocation and no locking, and reach the shared
 * per-reference accumulator in a single merge.
 *
 * PHMM_NO_PATH counts the read as rejected and contributes nothing else, the model having
 * no alignment of it to score. Anything else but PHMM_OK ends the run: the marginal then
 * fails on what would fail again on every read after it.
 *
 * target_pairs takes the read's co-modification, and is NULL for a run counting none. */
phmm_status tally(const cm_bam_record *read, const cm_fasta_record *ref,
                  const tally_tables *tables, tally_scratch *scratch,
                  accum *target, pairs *target_pairs);
