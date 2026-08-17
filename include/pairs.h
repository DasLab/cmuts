/* pairs.h -- co-modification of two reference positions.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include "phmm.h"

/* The sums one pair of positions carries, each a posterior expectation over reads.
 *
 * The first four give the 2x2 table a correlation is taken from: n11 = BOTH,
 * n10 = LEFT - BOTH, n01 = RIGHT - BOTH, and n00 = whatever of SPAN is left. The last is
 * not part of the table: the coverage is what the output reports as the reads behind a
 * pair. */
typedef enum {
    PAIR_SPAN,      /* both positions reached */
    PAIR_LEFT,      /* the lower modified, the upper reached */
    PAIR_RIGHT,     /* the upper modified, the lower reached */
    PAIR_BOTH,      /* both modified */
    PAIR_COVERED,   /* a base read at both, weighted by the quality of each */
    PAIR_N_CELLS,
} pair_cell;

/* The statistics derived from the pairs, as bits so a run may ask for any subset. What
 * is accumulated does not depend on the choice; only what is written does. */
typedef enum {
    PAIRS_CORRELATION = 1 << 0,   /* Pearson correlation of the two being modified */
    PAIRS_CONDITIONAL = 1 << 1,   /* probability one was modified where the other was */
} pairs_statistic;

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
 * The two positions are counted as independent given the read, which does not hold where
 * one path deviation spans both. */
void pairs_count(pairs *p, size_t len, const phmm_window *window);

/* Writes position i against every position: the correlation of the two being modified.
 * It is not a number where the reads are too few or either position is modified in all
 * of them or in none.
 *
 * The diagonal falls short of one wherever the posteriors are uncertain, and every
 * correlation involving that position is reduced the same way. */
void pairs_correlation(const pairs *p, size_t len, double min_depth, size_t i,
                       double *row);

/* Writes position i against every position: the probability i was modified in a read
 * where that position was, over the reads spanning both.
 *
 * The square is not symmetric: an entry conditions on the position it is indexed by. A
 * pair whose reads are too few, or whose conditioned position was never modified, is not
 * a number. */
void pairs_conditional(const pairs *p, size_t len, double min_depth, size_t i,
                       double *row);

void pairs_coverage(const pairs *p, size_t len, size_t i, double *row);
