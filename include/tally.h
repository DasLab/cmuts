/* tally.h -- one read's contribution to a reference's accumulated values.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "accum.h"
#include "bam.h"
#include "fasta.h"
#include "phmm.h"
#include "phred.h"

/* The settings a caller supplies, in types a command line can write into directly. */
typedef struct {
    int          band;     /* reference positions the marginal may look either side of the
                              CIGAR; 0 pins it to the path as written */
    phmm_weights weights;  /* what an event of each kind counts for */
} tally_config;

tally_config tally_defaults(void);

/* What the tally works from besides the read itself: the quality table and the model.
 *
 * Derived from the config rather than set alongside it, and built before any worker starts
 * and never written afterwards, so every worker may be given the same one. */
typedef struct {
    phred quality;
    phmm  model;
    int   band;   /* the half-width every row of the marginal is given */
} tally_tables;

void tally_tables_build(tally_tables *tables, const tally_config *config);

/* Working buffers one worker reuses across every read it processes.
 *
 * Private to a worker rather than shared, the marginal writing over the whole of it for
 * each read, and reused rather than allocated per read. Every read is counted through one,
 * so a caller must create it before any read arrives. */
typedef struct tally_scratch tally_scratch;

tally_scratch *tally_scratch_create(void);
void           tally_scratch_destroy(tally_scratch *scratch);

/* Adds one read's contribution to target, which is never cleared here.
 *
 * Accumulating into a caller-provided target rather than returning a fresh struct lets a
 * worker hold a private accumulator across a run of reads: contributions land there with
 * no allocation and no locking, and reach the shared per-reference accumulator in a single
 * merge.
 *
 * Anything but PHMM_OK ends the run rather than skipping the read, the marginal failing
 * only on what would fail again on every read after it. */
phmm_status tally(const cm_bam_record *read, const cm_fasta_record *ref,
                  const tally_tables *tables, tally_scratch *scratch,
                  accum *target);
