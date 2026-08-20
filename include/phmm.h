/* phmm.h -- marginalizing a read over the alignments a band admits.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "bam.h"
#include "fasta.h"
#include "params.h"
#include "phred.h"

/* How far either side of the CIGAR path the model may look. 2 admits every placement of a
 * one or two base deletion, and wider scored no better on the two libraries tested. */
#define PHMM_DEFAULT_BAND 2

/* The kinds of event the mutation channel counts. */
typedef enum {
    PHMM_SUBSTITUTION,
    PHMM_DELETION,
    PHMM_INSERTION,
    PHMM_N_EVENTS,
} phmm_event;

/* How much an event of each kind counts towards the mutation total, given that the event
 * happened. Read by the accumulation only, so the alignment does not depend on them. The
 * insertion weight also scales the span an insertion adds, so a weight of zero removes
 * insertions from both totals. */
typedef struct {
    double weight[PHMM_N_EVENTS];
} phmm_weights;

/* Substitutions and deletions at 1, insertions at 0. */
phmm_weights phmm_default_weights(void);

/* The parameters, with every transition they imply computed once. Insertion to deletion
 * and deletion to insertion are absent: an inserted base followed by a deleted one is a
 * substitution. A model is read-only once built, so every thread may share one. */
typedef struct {
    phmm_params  params;
    phmm_weights weights;
    double       match_to_match;
    double       match_to_insertion;
    double       match_to_deletion;
    double       insertion_to_insertion;
    double       insertion_to_match;
    double       deletion_to_deletion;
    double       deletion_to_match;
} phmm;

void phmm_build(phmm *model, const phmm_params *params,
                const phmm_weights *weights);

/* Buffers one thread reuses across reads, grown to fit the longest read seen. */
typedef struct phmm_scratch phmm_scratch;

phmm_scratch *phmm_scratch_create(void);
void          phmm_scratch_destroy(phmm_scratch *scratch);

/* What one read contributes to each reference position it could have reached.
 *
 * Values are fractional, spread over each placement the band allowed, and are added to a
 * target, never assigned. The window may extend past either end of the reference; those
 * positions are never written, and phmm_window_bounds gives the range inside it. */
typedef struct {
    hts_pos_t     origin;      /* reference position of value 0 */
    size_t        len;
    const double *coverage;    /* base read there, weighted by its quality */
    const double *spanned;     /* position reached, read or deleted */
    const double *mutations;   /* events attributed to that position */
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
