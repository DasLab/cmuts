/* pairs.h -- co-modification of two reference positions.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include "phmm.h"

/* The sums one pair of positions carries, each a posterior expectation over reads. They
 * give the 2x2 table a correlation is taken from: n11 = BOTH, n10 = LEFT - BOTH,
 * n01 = RIGHT - BOTH, and n00 = whatever of SPAN is left. */
typedef enum {
    PAIR_SPAN,    /* both positions reached */
    PAIR_LEFT,    /* the lower modified, the upper reached */
    PAIR_RIGHT,   /* the upper modified, the lower reached */
    PAIR_BOTH,    /* both modified */
    PAIR_N_CELLS,
} pair_cell;

/* Co-modification accumulated for one reference.
 *
 * Only i <= j is held, a pair being unordered. The packing is over the reference's own
 * length and not the capacity, so two may be merged only at the length both were zeroed
 * at. */
typedef struct {
    double *cells;   /* owned; PAIR_N_CELLS per pair, and NULL where none are held */
} pairs;

int  pairs_alloc(pairs *p, size_t cap);
void pairs_free(pairs *p);
void pairs_zero(pairs *p, size_t len);
void pairs_add(pairs *dst, const pairs *src, size_t len);

/* Adds one read's window, over the positions of it inside the reference.
 *
 * The two positions are taken as independent given the read, which holds beyond the reach
 * of a path deviation and not within it. */
void pairs_count(pairs *p, size_t len, const phmm_window *window);

/* Writes position i against every position: the correlation of the two being modified,
 * and the reads behind it.
 *
 * The coefficient is (BOTH * SPAN - LEFT * RIGHT) over the square root of the product of
 * the four marginals, and is not a number where the reads are too few or either position
 * is modified in all of them or in none.
 *
 * The diagonal falls short of one by however much the posteriors leave unsettled, which is
 * what every correlation involving that position is reduced by. */
void pairs_correlation(const pairs *p, size_t len, double min_depth, size_t i,
                       double *row);
void pairs_coverage(const pairs *p, size_t len, size_t i, double *row);
