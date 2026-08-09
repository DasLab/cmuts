/* phmm.h -- an alignment marginalized rather than taken at its word.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "bam.h"
#include "fasta.h"
#include "phred.h"

/* How far either side of the aligner's path the model may look.
 *
 * The band is drawn around the CIGAR and not around a diagonal, so what it
 * bounds is the departure from the alignment already found rather than the
 * departure from an ungapped one.
 *
 * How far it has to reach is set by the length of the gap and not by the length
 * of the run the gap sits in. Sliding a deletion of n bases along a homopolymer
 * puts every row it passes n positions off the path, however far it slides, so
 * two holds every placement of a deletion of one or two bases and a run of ten
 * asks for no more than a run of three. Wider admits only alignments that open
 * a gap twice or set a base against one it does not match, and on the two
 * libraries this was measured against -- which differ by nearly twofold in how
 * often they delete, and threefold in read length -- it scored no better and
 * mostly worse. Two is where both put the mass of their deletion lengths, and
 * a library whose gaps ran longer would want more.
 *
 * A row covers what the CIGAR path crosses on it before the band widens it at
 * all, so no band is ever too narrow to hold a read's own gaps and none is
 * refused for want of room. Nothing is the CIGAR path and nothing else, which
 * is the alignment as written and marginalized over alone.
 *
 * The maximum is what a command line will accept, the cost of a read being
 * linear in the band and nothing stopping a caller from asking for more than it
 * meant to. The model itself takes whatever it is given and asks for the memory
 * that implies. */
#define PHMM_DEFAULT_BAND 2
#define PHMM_MAX_BAND     256

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

/* The kinds of event the mutation channel counts. */
typedef enum {
    PHMM_SUBSTITUTION,
    PHMM_DELETION,
    PHMM_INSERTION,
    PHMM_N_EVENTS,
} phmm_event;

/* How much an event of each kind says about the template having been modified,
 * given that the event happened at all.
 *
 * Not a belief about how a read comes about, so nothing here reaches either
 * pass and the alignment does not depend on it: an insertion weighed at nothing
 * is still an insertion the read may have, and still the reason a spurious base
 * is not counted as a substitution instead. What a weight weighs is an event
 * already posited.
 *
 * Relative, in that a factor common to all of them rescales every rate the run
 * reports and changes nothing else. Whether any one of them is held at one is
 * the caller's to decide; they arrive at one apiece, which is the statement that
 * an event of any kind is worth one event.
 *
 * Measured from nothing, as the parameters are, and to be calibrated in time
 * against the enrichment each kind shows over an untreated control.
 *
 * Indexed by kind so that a finer weighing -- by reference base, say, the
 * chemistry not treating the four alike -- is another index rather than another
 * field. */
typedef struct {
    double weight[PHMM_N_EVENTS];
} phmm_weights;

phmm_weights phmm_default_weights(void);

/* The parameters with every transition they imply worked out once.
 *
 * A step from an insertion to a deletion is not among them, and neither is its
 * reverse: an inserted base followed by a deleted one is a substitution in
 * disguise, and allowing the pair would put the posterior for one event in two
 * places. Nothing writes to the model once built, so one may be shared by every
 * thread. */
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

/* How a marginalization ended.
 *
 * Neither failure is the read's doing, and neither leaves anything for a caller
 * to salvage. A row that sums to nothing takes a band with no width to it, a
 * row that sums to something not finite takes a parameter that is not a
 * probability, and the two passes disagreeing about which paths exist takes an
 * index gone wrong; none is a property of the alignment in hand, and none is
 * mended by trying the next one. They are told apart only so that a run ending
 * on one can say which it was. */
typedef enum {
    PHMM_OK,          /* out holds what the read is worth */
    PHMM_NO_MEMORY,   /* the matrix could not be had */
    PHMM_UNSOUND,     /* a pass did not hold together */
} phmm_status;

/* Marginalizes one read over the alignments its band admits, leaving the result
 * in out, which borrows from scratch and lasts until the next call on it.
 *
 * The band is given row by row: half[i] is how far either side of the CIGAR row
 * i may look, a row being one base of the placed span and one before them all.
 * It must hold at least read->l_qseq + 1 entries, that being the most rows a
 * read can have, and only those the span reaches are read. Constant throughout
 * is a band of one width; nothing here decides the shape.
 *
 * The read must store a sequence and place at least one of its bases, the
 * filter having turned away any that does neither.
 *
 * Anything but PHMM_OK ends the run: there is no read this can fail on and
 * leave the rest worth counting. Nothing is turned away for being large, either
 * -- what a read costs is what it costs, and it is the memory running out that
 * says so rather than a number chosen here. */
phmm_status phmm_run(const phmm *model, const phred *quality,
                     const cm_bam_record *read, const cm_fasta_record *ref,
                     const int *half, phmm_scratch *scratch, phmm_window *out);
