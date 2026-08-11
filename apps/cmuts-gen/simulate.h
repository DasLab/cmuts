/* simulate.h -- alignments built from a reference and a mutation model.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include <htslib/sam.h>

#include "sample.h"

/* How a read is drawn from a reference and made to differ from it.
 *
 * Insertions and deletions are counted as events per read with their own length rather than
 * falling out of a per-base rate. A single rate produces many one and two base indels and
 * effectively never a long one, where real aligners emit both -- and it is the long insertion,
 * leaving a read far longer than the reference it aligns to, that an upper length bound exists
 * to remove. */
typedef struct {
    distribution length;            /* reference span the read covers */
    double       mismatch_rate;     /* per aligned base */
    distribution insertions;        /* insertion events per read */
    distribution insertion_length;
    distribution deletions;
    distribution deletion_length;
    distribution soft_clips;        /* clipped ends, none through both */
    distribution soft_clip_length;
    distribution mapq;
    distribution base_quality;
    double       reverse_fraction;
} sim_model;

/* Buffers reused across reads, sized to the largest the model can produce. */
typedef struct sim_scratch sim_scratch;

sim_scratch *sim_scratch_create(const sim_model *model);
void         sim_scratch_destroy(sim_scratch *scratch);

/* Writes len random bases. */
void sim_sequence(char *out, size_t len, rng *r);

/* Where on a reference a read falls. Drawn separately from building the read so that a caller
 * can sort placements -- eight bytes apiece -- and emit the records in coordinate order, rather
 * than holding every record of a reference in memory to sort afterwards. */
typedef struct {
    hts_pos_t start;  /* 0-based */
    size_t    span;   /* reference bases covered */
} sim_placement;

sim_placement sim_place(const sim_model *model, rng *r, size_t reflen);

/* Builds one alignment at the given placement. The CIGAR, SEQ, MD and NM are derived from a
 * single list of per-base events, so they cannot disagree. Returns 0, or -1 if the record could
 * not be built. */
int sim_alignment(bam1_t *rec, sim_scratch *scratch, const sim_model *model,
                  rng *r, const char *name, int32_t tid,
                  const char *reference, sim_placement where);

/* Builds a record aligning nowhere, belonging to no reference. */
int sim_unmapped(bam1_t *rec, sim_scratch *scratch, const sim_model *model,
                 rng *r, const char *name);
