/* phmm.h -- marginalizing a read over the alignments a band admits.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "bam.h"
#include "fasta.h"
#include "phred.h"

/* How far either side of the aligner's path the model may look.
 *
 * The band is drawn around the CIGAR, not around a diagonal, so it bounds
 * departure from the reported alignment. Required width follows the length of
 * the gap, not the length of the homopolymer run the gap sits in: 2 admits
 * every placement of a one or two base deletion. Wider scored no better on the
 * two libraries tested.
 *
 * A row always covers what the CIGAR path crosses, before the band widens it,
 * so no band is too narrow to hold a read's own gaps. A band of 0 is the CIGAR
 * path alone.
 *
 * There is no upper bound. Cost per read is linear in the band, and an
 * unreasonable one ends the run by exhausting memory. */
#define PHMM_DEFAULT_BAND 2

/* Provisional default rates, the right order of magnitude for reverse
 * transcription of a chemically modified template but not measured. They should
 * become command-line options, and eventually estimates from the alignments. */
typedef struct {
    double open_deletion;     /* chance an aligned position begins one */
    double open_insertion;
    double extend_deletion;   /* geometric, so 1 - 1 / mean length */
    double extend_insertion;
    double modification;      /* chance a base differs from the reference for
                                 some reason other than a misread */
} phmm_params;

phmm_params phmm_defaults(void);

/* The kinds of event the mutation channel counts. */
typedef enum {
    PHMM_SUBSTITUTION,
    PHMM_DELETION,
    PHMM_INSERTION,
    PHMM_N_EVENTS,
} phmm_event;

/* How much an event of each kind counts towards the mutation total, given that
 * the event happened.
 *
 * These reach the accumulation only. Neither pass reads them, so the alignment
 * does not depend on them. They are relative: a factor common to all three
 * rescales every rate the run reports and changes nothing else. Each defaults
 * to 1, and they are uncalibrated.
 *
 * Indexed by kind so that a finer breakdown, by reference base for example, is
 * another index rather than another field. */
typedef struct {
    double weight[PHMM_N_EVENTS];
} phmm_weights;

phmm_weights phmm_default_weights(void);

/* The parameters with every transition they imply worked out once.
 *
 * Insertion to deletion and deletion to insertion are absent: an inserted base
 * followed by a deleted one is a substitution, and allowing the pair would
 * record one event's posterior in two places. The model is read-only once
 * built, so one may be shared by every thread. */
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
 * A window, not the whole reference, since a read reaches a few hundred
 * positions of a reference that may be far longer. Values are fractional,
 * spread over each placement the band allowed, and must be added to a target
 * rather than assigned.
 *
 * The window is sized from the band, which is not clamped, so it may extend
 * past either end of the reference. Nothing is ever laid there, but the
 * positions are present, and the caller must keep them out of anything it
 * indexes by reference position. */
typedef struct {
    hts_pos_t     origin;      /* reference position of value 0 */
    size_t        len;
    const double *coverage;    /* base read there, weighed by its quality */
    const double *spanned;     /* position reached, read or deleted */
    const double *mutations;   /* events attributed to that position */
} phmm_window;

/* How a marginalization ended.
 *
 * Both failures indicate a bad parameter or a bug, not a bad read, and neither
 * is fixed by trying the next alignment. A row summing to zero means a band
 * with no width; a row summing to something non-finite means a parameter that
 * is not a probability; the two passes disagreeing about which paths exist
 * means an index error. They are distinguished only so a run can report which
 * occurred. */
typedef enum {
    PHMM_OK,          /* out holds the read's contribution */
    PHMM_NO_MEMORY,   /* the matrix could not be allocated */
    PHMM_UNSOUND,     /* a pass did not hold together */
} phmm_status;

/* Marginalizes one read over the alignments its band admits, leaving the result
 * in out, which borrows from scratch and is valid until the next call on it.
 *
 * The band is given per row: half[i] is how far either side of CIGAR row i the
 * model may look, a row being one base of the placed span plus one before them
 * all. half must hold at least read->l_qseq + 1 entries; only those the span
 * reaches are read. A constant array gives a uniform band. The shape is the
 * caller's to choose.
 *
 * The read must store a sequence and place at least one base; the filter
 * rejects any that does neither.
 *
 * Anything but PHMM_OK ends the run. No read is rejected for being large: cost
 * is whatever the read implies, and memory exhaustion is what reports it. */
phmm_status phmm_run(const phmm *model, const phred *quality,
                     const cm_bam_record *read, const cm_fasta_record *ref,
                     const int *half, phmm_scratch *scratch, phmm_window *out);
