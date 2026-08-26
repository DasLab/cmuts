/* phmm.h -- a banded pair HMM.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "bam.h"
#include "fasta.h"
#include "params.h"
#include "phred.h"

/* How far either side of the CIGAR path the model may look. */
#define PHMM_DEFAULT_BAND 2

/* The transition weights that hold across every reference: the extending of each run
 * and the closing that leaves what the extension does not take. The modification and
 * opening rates vary per base and arrive through the profile instead. Insertion to
 * deletion and deletion to insertion are absent: an inserted base followed by a
 * deleted one is a substitution. */
typedef struct {
    double match_to_insertion;
    double match_to_deletion;
    double insertion_to_insertion;
    double deletion_to_deletion;
} phmm;

void phmm_build(phmm *model, const phmm_params *params);

/* Per-base rates, one value per reference base, indexed as the counted events are.
 * modification[b] is the chance base b was modified before it was read.
 * open_insertion[b] opens the insertion counted at base b, and open_deletion[b] the
 * deletion run whose 3'-most deleted base is b; the two are the arms of the decision
 * made at the pairing 3' of each run. A cell outside the reference reads the nearest
 * value, which reaches only terms a zero emission clears. */
typedef struct {
    const double *modification;
    const double *open_insertion;
    const double *open_deletion;
} phmm_profile;

/* Buffers one thread reuses across reads, grown to fit the longest read seen. */
typedef struct phmm_scratch phmm_scratch;

phmm_scratch *phmm_scratch_create(void);
void          phmm_scratch_destroy(phmm_scratch *scratch);

/* Computed values for one read, in a window around the reference. */
typedef struct {
    hts_pos_t     origin;       /* reference position of value 0 */
    size_t        len;
    const double *coverage;     /* base read there */
    const double *mismatches;   /* template differences under the read bases */
    const double *insertions;   /* insertions opened after this base */
    const double *deletions;    /* deletions opened at this base */
    const double *ends;         /* the read's 5'-most pairing lands here */
} phmm_window;

/* Gives the window indices that fall inside a reference of len bases: from the first to
 * one past the last, an empty range where none do. */
void phmm_window_bounds(const phmm_window *window, size_t len,
                        size_t *from, size_t *to);

/* How a marginalization ended. PHMM_NO_PATH skips the read, while the other
 * two failures end the run. */
typedef enum {
    PHMM_OK,          /* out holds the read's contribution */
    PHMM_NO_MEMORY,   /* the matrix could not be allocated */
    PHMM_NO_PATH,     /* no alignment the band admits has any probability */
    PHMM_UNSOUND,     /* the two passes did not agree */
} phmm_status;

/* Marginalizes one read over the alignments its band admits, leaving the result
 * in out, which borrows from scratch and is valid until the next call on it.
 *
 * The band is given per row: half[i] is how far either side of CIGAR row i the model may
 * look, a row being one base of the placed span plus one before them all. half must hold
 * at least read->l_qseq + 1 entries, of which only those the span reaches are read.
 *
 * Each profile array must hold one value per base of ref.
 *
 * The read must store a sequence and place at least one base, which the filter
 * enforces. */
phmm_status phmm_run(const phmm *model, const phmm_profile *profile,
                     const phred *quality, const cm_bam_record *read,
                     const cm_fasta_record *ref, const int *half,
                     phmm_scratch *scratch, phmm_window *out);
