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

/* Transition rates. Insertion to deletion and deletion to insertion are absent:
 * an inserted base followed by a deleted one is a substitution. */
typedef struct {
    phmm_params  params;
    double       match_to_match;
    double       match_to_insertion;
    double       match_to_deletion;
    double       insertion_to_insertion;
    double       insertion_to_match;
    double       deletion_to_deletion;
    double       deletion_to_match;
} phmm;

void phmm_build(phmm *model, const phmm_params *params);

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
 * The read must store a sequence and place at least one base, which the filter
 * enforces. */
phmm_status phmm_run(const phmm *model, const phred *quality,
                     const cm_bam_record *read, const cm_fasta_record *ref,
                     const int *half, phmm_scratch *scratch, phmm_window *out);
