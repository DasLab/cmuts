/* refrow.c -- each output field's values, taken from an accumulator.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "refrow.h"

#include <stdlib.h>

#include "output.h"

struct refrow {
    h5writer   *out;      /* borrowed */
    rate_config rates;
    double     *row;      /* scratch a derived field is computed into */
};

refrow *refrow_create(h5writer *out, rate_config rates, size_t ref_cap,
                      const bool *wanted)
{
    refrow *r      = calloc(1, sizeof *r);
    size_t  widest = out_widest(ref_cap, wanted);

    if (!r) {
        return NULL;
    }

    r->out   = out;
    r->rates = rates;
    r->row   = calloc(widest ? widest : 1, sizeof *r->row);

    if (!r->row) {
        refrow_destroy(r);
        return NULL;
    }

    return r;
}

void refrow_destroy(refrow *r)
{
    if (!r) {
        return;
    }

    free(r->row);
    free(r);
}

/* Fills the scratch with one value for every ordered pair, a row of the reference at a
 * time. What is accumulated is the triangle; the square is written whole, so a reader
 * indexes it either way round. Only the conditional differs with the order, an entry
 * conditioning on the position of its column. */
static void pair_square(refrow *r, out_field_id id, const pairs *pr, size_t len)
{
    for (size_t i = 0; i < len; i++) {
        double *row = r->row + i * len;

        if (id == OUT_PAIRWISE_CORRELATION) {
            pairs_correlation(pr, len, r->rates.min_depth, i, row);
        } else if (id == OUT_PAIRWISE_CONDITIONAL) {
            pairs_conditional(pr, len, r->rates.min_depth, i, row);
        } else {
            pairs_coverage(pr, len, i, row);
        }
    }
}

/* Gives one output field's values, computed into the scratch row where they are derived
 * and read in place where they are not, or NULL for a field with no row.
 * The accumulated fields and the written ones do not correspond one to one, so every
 * field is listed here and none is defaulted: one added without a source of its own
 * draws a warning and is refused at the write. */
static const double *values(refrow *r, out_field_id id, const accum *acc,
                            const pairs *pr, size_t len)
{
    switch (id) {
        case OUT_COVERAGE:   return accum_const_data(acc, ACCUM_COVERAGE);
        case OUT_LENGTHS:    return accum_const_data(acc, ACCUM_LENGTHS);
        case OUT_READS:      return accum_const_data(acc, ACCUM_READS);
        case OUT_REJECTED:   return accum_const_data(acc, ACCUM_FILTERED);
        case OUT_NORM:
        case OUT_UNMAPPED:
        case OUT_N_FIELDS:   break;
        case OUT_REACTIVITY:
            rate_reactivity(&r->rates, acc, len, r->row);
            return r->row;
        case OUT_ERROR:
            rate_error(&r->rates, acc, len, r->row);
            return r->row;
        case OUT_PAIRWISE_CORRELATION:
        case OUT_PAIRWISE_CONDITIONAL:
        case OUT_PAIRWISE_COVERAGE:
            if (!pr) {
                break;
            }
            pair_square(r, id, pr, len);
            return r->row;
    }

    return NULL;
}

/* Whether a field's row spans more than one extent, and so is written as a block rather
 * than as a row with a tail to mark. */
static bool is_block(out_field_id id, size_t len)
{
    return shape_rank(OUT_FIELDS[id].row(len, len)) > 1;
}

int refrow_write(refrow *r, int32_t tid, size_t len, const accum *acc,
                 const pairs *pr)
{
    if (len == 0) {
        return 0;
    }

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        const double *row;

        if (!OUT_FIELDS[id].per_ref || !h5writer_holds(r->out, id)) {
            continue;
        }

        row = values(r, id, acc, pr, len);

        if (!row) {
            return -1;
        }

        if (is_block(id, len)
            ? h5writer_block(r->out, id, tid, len, row) < 0
            : h5writer_field(r->out, id, tid, len, row) < 0) {
            return -1;
        }
    }

    return 0;
}
