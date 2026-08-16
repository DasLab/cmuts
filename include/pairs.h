/* pairs.h -- co-modification of two reference positions.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include "phmm.h"

/* The sums one pair of positions carries.
 *
 * Each is a posterior expectation summed over reads, and together they give the table a
 * correlation is taken from: n11 = BOTH, n10 = LEFT - BOTH, n01 = RIGHT - BOTH, and
 * n00 = SPAN - LEFT - RIGHT + BOTH, which sum to SPAN. */
typedef enum {
    PAIR_SPAN,    /* both positions reached */
    PAIR_LEFT,    /* the lower modified, the upper reached */
    PAIR_RIGHT,   /* the upper modified, the lower reached */
    PAIR_BOTH,    /* both modified */
    PAIR_N_CELLS,
} pair_cell;

/* Co-modification accumulated for one reference.
 *
 * A pair is unordered, so only i <= j is held. The packing is over the reference's own
 * length and not the capacity every accumulator is allocated at, so two of these may be
 * merged only at the length both were zeroed at. Every caller already works that way: a
 * shadow and the reference it merges into are zeroed at the same length. */
typedef struct {
    double *cells;   /* owned; PAIR_N_CELLS per pair */
    size_t  cap;
} pairs;

int  pairs_alloc(pairs *p, size_t cap);
void pairs_free(pairs *p);
void pairs_zero(pairs *p, size_t len);
void pairs_add(pairs *dst, const pairs *src, size_t len);

/* Adds one read's window, over the positions of it that fall inside the reference.
 *
 * The two positions are taken as independent given the read, which they are beyond the
 * reach of a path deviation and are not within it. */
void pairs_count(pairs *p, size_t len, const phmm_window *window);

/* Writes position i against every position of the reference: the correlation of the two
 * being modified, and the reads behind it.
 *
 * The correlation is the Pearson coefficient of two binary variables, which for a table
 * of these four sums is (BOTH * SPAN - LEFT * RIGHT) over the square root of the product
 * of the four marginals. It is not a number where the reads are too few, where either
 * position is modified in all of them or in none, and on the diagonal, none of which
 * carry a correlation. */
void pairs_correlation(const pairs *p, size_t len, double min_depth, size_t i,
                       double *row);
void pairs_coverage(const pairs *p, size_t len, size_t i, double *row);
