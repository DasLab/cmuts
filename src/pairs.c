/* pairs.c -- co-modification of two reference positions, accumulated per read.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "pairs.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* Pairs a reference of len bases holds, which is its upper triangle including the
 * diagonal. */
static size_t pair_count(size_t len)
{
    return len * (len + 1) / 2;
}

static size_t pair_values(size_t len)
{
    return pair_count(len) * PAIR_N_CELLS;
}

/* Offset of the pair (i, j), for i <= j, in a reference of len bases. Row i holds the
 * len - i pairs from the diagonal rightwards, so the rows before it hold i * (2 * len -
 * i + 1) / 2 pairs between them. */
static size_t pair_at(size_t len, size_t i, size_t j)
{
    return (i * (2 * len - i + 1) / 2 + (j - i)) * PAIR_N_CELLS;
}

static double *cell(pairs *p, size_t len, size_t i, size_t j)
{
    return p->cells + pair_at(len, i, j);
}

static const double *const_cell(const pairs *p, size_t len, size_t i, size_t j)
{
    return p->cells + pair_at(len, i, j);
}

/* ------------------------------------------------------------------------ */
/* Storage                                                                   */
/* ------------------------------------------------------------------------ */

int pairs_alloc(pairs *p, size_t cap)
{
    size_t values = pair_values(cap);

    p->cells = calloc(values ? values : 1, sizeof *p->cells);

    return p->cells ? 0 : -1;
}

void pairs_free(pairs *p)
{
    free(p->cells);
    p->cells = NULL;
}

void pairs_zero(pairs *p, size_t len)
{
    memset(p->cells, 0, pair_values(len) * sizeof *p->cells);
}

void pairs_add(pairs *dst, const pairs *src, size_t len)
{
    size_t n = pair_values(len);

    for (size_t i = 0; i < n; i++) {
        dst->cells[i] += src->cells[i];
    }
}

/* ------------------------------------------------------------------------ */
/* Counting                                                                  */
/* ------------------------------------------------------------------------ */

void pairs_count(pairs *p, size_t len, const phmm_window *window)
{
    size_t from;
    size_t to;

    phmm_window_bounds(window, len, &from, &to);

    for (size_t a = from; a < to; a++) {
        size_t i  = (size_t)(window->origin + (hts_pos_t)a);
        double si = window->spanned[a];
        double mi = window->mutations[a];
        double ci = window->coverage[a];

        for (size_t b = a; b < to; b++) {
            size_t  j  = (size_t)(window->origin + (hts_pos_t)b);
            double  sj = window->spanned[b];
            double  mj = window->mutations[b];
            double *at = cell(p, len, i, j);

            at[PAIR_SPAN]    += si * sj;
            at[PAIR_LEFT]    += mi * sj;
            at[PAIR_RIGHT]   += si * mj;
            at[PAIR_BOTH]    += mi * mj;
            at[PAIR_COVERED] += ci * window->coverage[b];
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Reading                                                                   */
/* ------------------------------------------------------------------------ */

/* Gives the pair (i, j) in either order. */
static const double *ordered(const pairs *p, size_t len, size_t i, size_t j)
{
    return i <= j ? const_cell(p, len, i, j) : const_cell(p, len, j, i);
}

/* Gives the correlation of the four sums, or NaN where they support none: the marginals
 * are the modified and unmodified reads at each position, and a coefficient needs all
 * four positive. The formula is symmetric in the two positions, so the pair may be
 * fetched in either order. */
static double coefficient(const double *at, double min_depth)
{
    double span  = at[PAIR_SPAN];
    double left  = at[PAIR_LEFT];
    double right = at[PAIR_RIGHT];
    double both  = at[PAIR_BOTH];
    double scale = left * (span - left) * right * (span - right);

    if (span <= 0.0 || span < min_depth || scale <= 0.0) {
        return (double)NAN;
    }

    return (both * span - left * right) / sqrt(scale);
}

void pairs_correlation(const pairs *p, size_t len, double min_depth, size_t i,
                       double *row)
{
    for (size_t j = 0; j < len; j++) {
        row[j] = coefficient(ordered(p, len, i, j), min_depth);
    }
}

/* Gives the probability the row's position was modified in a read where the conditioned
 * position was, or NaN where the sums support none. */
static double conditional(const double *at, double conditioned, double min_depth)
{
    double span = at[PAIR_SPAN];

    if (span <= 0.0 || span < min_depth || conditioned <= 0.0) {
        return (double)NAN;
    }

    return at[PAIR_BOTH] / conditioned;
}

void pairs_conditional(const pairs *p, size_t len, double min_depth, size_t i,
                       double *row)
{
    for (size_t j = 0; j < len; j++) {
        const double *at = ordered(p, len, i, j);

        /* For i <= j the stored pair holds j as its upper, whose mutations are
         * PAIR_RIGHT; otherwise j is the lower and they are PAIR_LEFT. */
        row[j] = conditional(at, i <= j ? at[PAIR_RIGHT] : at[PAIR_LEFT], min_depth);
    }
}

void pairs_coverage(const pairs *p, size_t len, size_t i, double *row)
{
    for (size_t j = 0; j < len; j++) {
        row[j] = ordered(p, len, i, j)[PAIR_COVERED];
    }
}
