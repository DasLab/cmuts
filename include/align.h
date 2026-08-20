/* align.h -- where a CIGAR puts each base of a read.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdint.h>

#include "bam.h"

/* The stretch of a read the alignment places on the reference. Clipped bases are placed
 * nowhere and lie outside it. */
typedef struct {
    int32_t begin;  /* first query offset the alignment places */
    int32_t end;    /* one past the last */
} aln_span;

/* Gives the placed span, empty where the alignment places no base: a record carrying no
 * CIGAR, one whose every operation is a clip, or one consuming reference but no read. */
aln_span aln_placed_span(const cm_bam_record *read);

/* How far into the reference the CIGAR has reached at one point of the placed span: first
 * once that base has been taken, last once the deletion straight after it has been too.
 * The two differ only where a deletion follows. */
typedef struct {
    hts_pos_t first;
    hts_pos_t last;
} aln_place;

/* Fills places with one for each point of the placed span, and returns the span itself:
 * places[0] is where the read starts, and one follows for every base it places. Places
 * are indexed by prefix length. places must hold read->l_qseq + 1 values. */
aln_span aln_places(const cm_bam_record *read, aln_place *places);
