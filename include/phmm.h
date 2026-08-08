/* phmm.h -- an alignment marginalized rather than taken at its word.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>

#include "bam.h"
#include "fasta.h"
#include "phred.h"

/* How far either side of the aligner's path the model may look.
 *
 * The band is drawn around the CIGAR and not around a diagonal, so what it
 * bounds is the departure from the alignment already found rather than the
 * departure from an ungapped one. Sixteen is far wider than the homopolymer
 * runs that make a deletion's position ambiguous in the first place, and the
 * cost of a read is linear in it. */
#define PHMM_BAND 16

/* Reads longer than this are left to the walk. The matrix a read needs is its
 * length times the band times the three states, and a read long enough for that
 * to matter is long enough that its alignment is a different problem. */
#define PHMM_MAX_QUERY 20000

/* What the model believes before it has seen a read.
 *
 * Provisional. These are the right order of magnitude for reverse transcription
 * of a chemically modified template and are measured from nothing; they belong
 * on the command line, and ought in time to be estimated from the alignments
 * themselves. Until then they sit here, in one place, so that there is a single
 * thing to change. */
typedef struct {
    double open_deletion;     /* chance an aligned position begins one */
    double open_insertion;
    double extend_deletion;   /* geometric, so 1 - 1 / mean length */
    double extend_insertion;
    double modification;      /* chance a base differs from the reference for
                                 some reason other than a misread */
} phmm_params;

phmm_params phmm_defaults(void);

/* The parameters with every transition they imply worked out once.
 *
 * A step from an insertion to a deletion is not among them, and neither is its
 * reverse: an inserted base followed by a deleted one is a substitution in
 * disguise, and allowing the pair would put the posterior for one event in two
 * places. Nothing writes to the model once built, so one may be shared by every
 * thread. */
typedef struct {
    phmm_params params;
    double      match_to_match;
    double      match_to_insertion;
    double      match_to_deletion;
    double      insertion_to_insertion;
    double      insertion_to_match;
    double      deletion_to_deletion;
    double      deletion_to_match;
} phmm;

void phmm_build(phmm *model, const phmm_params *params);

/* Having read what was read, the chance the template really differed from the
 * reference here.
 *
 * A base agreeing with the reference may have been modified and then misread
 * back into agreement, and one disagreeing may be an unmodified base misread.
 * The first is negligible and the second is not, which is what makes a poorly
 * read disagreement worth less than a clean one.
 *
 * Exposed because the walk counts the same quantity the marginal does, and the
 * two have to agree wherever the alignment was never in doubt. */
double phmm_modification(const phmm *model, bool agree, double error);

/* Buffers one thread reuses across reads, grown to whatever the longest read so
 * far has needed. */
typedef struct phmm_scratch phmm_scratch;

phmm_scratch *phmm_scratch_create(void);
void          phmm_scratch_destroy(phmm_scratch *scratch);

/* What one read is worth to each reference position it could have reached.
 *
 * A window rather than the whole reference, because a read reaches a few
 * hundred positions of a reference that may be far longer. The values are the
 * read's alone and are to be added, not assigned, and every one of them is
 * fractional, being spread over each placement the band allowed.
 *
 * The window may begin before the reference does or end after it: a read placed
 * near either boundary really has paths that leave it, and dropping them is the
 * caller's to do. */
typedef struct {
    hts_pos_t     origin;      /* reference position of value 0 */
    size_t        len;
    const double *coverage;    /* base read there, weighed by its quality */
    const double *spanned;     /* position reached, read or deleted */
    const double *mutations;   /* events laid at that position's door */
} phmm_window;

/* Marginalizes one read over the alignments its band admits, leaving the result
 * in out, which borrows from scratch and lasts until the next call on it.
 *
 * Returns false where the read cannot be marginalized -- no sequence, nothing
 * placed, a read or a span too long, or a pass that came to nothing -- and the
 * caller is to fall back on the alignment as it was written. */
bool phmm_run(const phmm *model, const phred *quality,
              const cm_bam_record *read, const cm_fasta_record *ref,
              phmm_scratch *scratch, phmm_window *out);
