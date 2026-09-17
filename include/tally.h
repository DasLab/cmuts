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
    int                band;       /* reference positions the marginal may look either
                                      side of the CIGAR; 0 pins it to the path as
                                      written */
    int                min_phred;  /* a base scoring below this takes the maximum error */
    phmm_uniform_rates uniform;    /* the rates that fill every base of every reference */
} tally_config;

tally_config tally_defaults(void);

/* What the tally works from besides the read and its reference's rates: the quality
 * table, the rates every reference is filled from, and the band. Built once from the
 * config and never written afterwards, so every worker may share one. */
typedef struct {
    phred              quality;
    phmm_uniform_rates uniform;  /* the rates tally_rates_fill writes */
    int                band;     /* the half-width every row of the marginal is given */
} tally_tables;

void tally_tables_build(tally_tables *tables, const tally_config *config);

/* Fills the rates of a reference of len bases from the uniform rates, with ends equal at
 * every base. The arrays of rates borrow storage, which must hold PHMM_RATE_ARRAYS times
 * len values and outlive every use of rates. */
void tally_rates_fill(phmm_rates *rates, double *storage, size_t len,
                      const tally_tables *tables);

/* Working buffers one worker reuses across every read it processes. Private to that
 * worker; the marginal writes over the whole of it for each read. */
typedef struct tally_scratch tally_scratch;

tally_scratch *tally_scratch_create(void);
void           tally_scratch_destroy(tally_scratch *scratch);

/* Adds one read's contribution to target, which is never cleared here. PHMM_NO_PATH
 * adds one to the rejected count and no other value; any other failure ends the run.
 * target_pairs takes the read's co-modification, and is NULL for a run counting none.
 * model is prepared for ref with a margin of at least the band. */
phmm_status tally(const cm_bam_record *read, const cm_fasta_record *ref,
                  const phmm_model *model, const tally_tables *tables,
                  tally_scratch *scratch, accum *target, pairs *target_pairs);
