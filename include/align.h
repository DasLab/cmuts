/* align.h -- where a CIGAR puts each base of a read.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdint.h>

#include "bam.h"

/* The stretch of a read the alignment places on the reference.
 *
 * Clipped bases are placed nowhere, so anything free to move a base must be told which bases
 * it may move. Left unconstrained it would align the clipped ends too, inventing a placement
 * for bases the aligner left unplaced. */
typedef struct {
    int32_t begin;  /* first query offset the alignment places */
    int32_t end;    /* one past the last */
} aln_span;

/* The placed span, empty where the alignment places nothing: a record carrying no CIGAR, one
 * whose every operation is a clip, or one naming reference alone. Costs a look at either
 * end, both stopping at the first operation that is not a clip. */
aln_span aln_placed_span(const cm_bam_record *read);

/* How far into the reference the CIGAR has reached at one point of the placed span: first
 * once that base has been taken, last once whatever is skipped straight after has been too.
 *
 * The two are equal unless a deletion follows, and what lies between them is reference the
 * read reached without reading. A point is therefore a stretch and not a position, which is
 * the whole of what a deletion does to the path: it makes one base of the read account for
 * several of the reference. */
typedef struct {
    hts_pos_t first;
    hts_pos_t last;
} aln_place;

/* Fills places with one for each point of the placed span, and returns the span itself:
 * places[0] is where the read starts, and one follows for every base it places.
 *
 * Prefix lengths rather than coordinates, that being the indexing a dynamic program over two
 * sequences runs on, and leaving a deletion as a step along one of them rather than a hole in
 * the other. places must hold read->l_qseq + 1 values, which is always enough. */
aln_span aln_places(const cm_bam_record *read, aln_place *places);
