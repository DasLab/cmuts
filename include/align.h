/* align.h -- an alignment read as a sequence of runs.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stdint.h>

#include "bam.h"
#include "fasta.h"

/* What a stretch of one CIGAR operation turned out to be.
 *
 * An aligned operation says only that the read and the reference were placed
 * against each other, not that they agree -- and aligners write M rather than
 * = or X almost always -- so it is settled by comparing the two. A stretch of
 * M becomes runs of MATCH, MISMATCH and AMBIGUOUS.
 *
 * The others are what they say they are. Ambiguity belongs to a comparison and
 * not to a base: an inserted N is an insertion, and a deleted N is a deletion,
 * whatever letter either happens to be. */
typedef enum {
    ALN_MATCH,
    ALN_MISMATCH,
    ALN_AMBIGUOUS,   /* aligned, but one side is not a named base */
    ALN_INSERTION,
    ALN_DELETION,
    ALN_SKIP,        /* reference passed over, CIGAR N */
    ALN_SOFT_CLIP,
    ALN_HARD_CLIP,
    ALN_PAD,
} aln_kind;

/* Runs rather than single positions, so that what counts as one event is the
 * caller's to decide: three mismatches in a row arrive as one run of three. */
typedef struct {
    aln_kind  kind;
    hts_pos_t reference;  /* first reference position the run covers */
    int32_t   query;      /* first offset into the read */
    uint32_t  len;
} aln_run;

/* Every operation is reported, including those advancing neither sequence, so
 * that a caller reading only the runs stays in step with both. */
typedef struct {
    const cm_bam_record   *read;
    const cm_fasta_record *ref;
    uint32_t               op;        /* operation being reported */
    uint32_t               consumed;  /* how much of it already has been */
    hts_pos_t              reference;
    int32_t                query;
} aln_walk;

/* The read must store a sequence, the filter having turned away any that does
 * not: there is no comparing what is not there. */
void aln_open(aln_walk *walk, const cm_bam_record *read, const cm_fasta_record *ref);

/* Fills run with the next one, or returns false at the end of the read. */
bool aln_next(aln_walk *walk, aln_run *run);
