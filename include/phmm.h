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

/* Every parameter of the model for one reference.
 *
 * The transition weights hold along the whole reference: the extending of each run and
 * the closing that leaves what the extension does not take. Insertion to deletion and
 * deletion to insertion are absent: an inserted base followed by a deleted one is a
 * substitution.
 *
 * The arrays hold one value per reference base, indexed as the counted events are.
 * modification[b] is the chance base b was modified before it was read.
 * open_insertion[b] opens the insertion counted at base b, and open_deletion[b] the
 * deletion run whose 3'-most deleted base is b; the two are the arms of the decision
 * made at the pairing 3' of each run. ends[b] is proportional to the number of reads
 * whose 5'-most paired base is b, which is the termination rate at b times the coverage
 * at b; only the ratios between bases in one read's band matter. A cell outside the
 * reference reads the nearest value, which reaches only terms a zero emission clears. */
typedef struct {
    double        match_to_insertion;
    double        match_to_deletion;
    double        insertion_to_insertion;
    double        deletion_to_deletion;
    const double *modification;
    const double *open_insertion;
    const double *open_deletion;
    const double *ends;
} phmm_rates;

/* The arrays of phmm_rates. */
#define PHMM_RATE_ARRAYS 4

/* Sets the transition weights from the extension rates in uniform. The arrays are left
 * as they are. */
void phmm_rates_set_transitions(phmm_rates *rates, const phmm_uniform_rates *uniform);

/* One reference base as every cell that pairs it or steps into its pairing reads it.
 * Defined with the model. */
typedef struct phmm_column phmm_column;

/* One reference prepared for the marginal: its rates, and one column for each base
 * within margin bases of either end of it, so that a band reaching past an end still
 * finds a column. Built once per reference and read by every read against it, from any
 * thread. */
typedef struct {
    const phmm_rates *rates;
    phmm_column      *columns;  /* owned */
    size_t            cap;      /* columns allocated */
    size_t            len;      /* bases of the reference prepared */
    int               margin;   /* bases covered either side of the reference */
} phmm_model;

/* Prepares ref under rates, whose arrays must hold one value per base of ref and must
 * outlive every use of the model. The columns are grown to fit and never shrink, so a
 * model can be prepared again for another reference. Returns -1 where they cannot be
 * grown. */
int phmm_model_prepare(phmm_model *model, const phmm_rates *rates,
                       const cm_fasta_record *ref, int margin);

void phmm_model_free(phmm_model *model);

/* Buffers one thread reuses across reads, grown to fit the longest read seen. */
typedef struct phmm_scratch phmm_scratch;

phmm_scratch *phmm_scratch_create(void);
void          phmm_scratch_destroy(phmm_scratch *scratch);

/* The computed values for one read at one reference position. */
typedef struct {
    double coverage;     /* base read there */
    double mismatches;   /* template differences under the read bases */
    double insertions;   /* insertions opened after this base */
    double deletions;    /* deletions opened at this base */
    double ends;         /* the read's 5'-most pairing lands here */
} phmm_position;

/* Computed values for one read, in a window around the reference. */
typedef struct {
    hts_pos_t            origin;  /* reference position of value 0 */
    size_t               len;
    const phmm_position *at;      /* one per position of the window */
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
 * at least read->l_qseq + 1 entries, of which only those the span reaches are read, and
 * none may exceed the margin the model was prepared with.
 *
 * The model must be prepared for ref. A read placed beyond the bases the model covers
 * cannot be scored against it and is given PHMM_NO_PATH; an alignment within its
 * reference never is.
 *
 * The read must store a sequence and place at least one base, which the filter
 * enforces. */
phmm_status phmm_run(const phmm_model *model, const phred *quality,
                     const cm_bam_record *read, const cm_fasta_record *ref,
                     const int *half, phmm_scratch *scratch, phmm_window *out);
