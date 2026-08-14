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

refrow *refrow_create(h5writer *out, rate_config rates, size_t ref_cap)
{
    refrow *r = calloc(1, sizeof *r);

    if (!r) {
        return NULL;
    }

    r->out   = out;
    r->rates = rates;
    r->row   = calloc(ref_cap ? ref_cap : 1, sizeof *r->row);

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

/* Gives one output field's values, computed into the scratch row where they are derived
 * and read in place where they are not, or NULL for a field with no row.
 * The accumulated fields and the written ones do not correspond one to one, so every
 * field is listed here and none is defaulted: one added without a source of its own
 * draws a warning and is refused at the write. */
static const double *values(refrow *r, out_field_id id, const accum *acc,
                            size_t len)
{
    switch (id) {
        case OUT_COVERAGE:   return accum_const_data(acc, ACCUM_COVERAGE);
        case OUT_LENGTHS:    return accum_const_data(acc, ACCUM_LENGTHS);
        case OUT_READS:      return accum_const_data(acc, ACCUM_READS);
        case OUT_REJECTED:   return accum_const_data(acc, ACCUM_FILTERED);
        case OUT_UNMAPPED:
        case OUT_N_FIELDS:   break;
        case OUT_REACTIVITY:
            rate_reactivity(&r->rates, acc, len, r->row);
            return r->row;
        case OUT_ERROR:
            rate_error(&r->rates, acc, len, r->row);
            return r->row;
    }

    return NULL;
}

int refrow_write(refrow *r, int32_t tid, size_t len, const accum *acc)
{
    if (len == 0) {
        return 0;
    }

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        const double *row;

        if (!OUT_FIELDS[id].per_ref) {
            continue;
        }

        row = values(r, id, acc, len);

        if (!row || h5writer_field(r->out, id, tid, len, row) < 0) {
            return -1;
        }
    }

    return 0;
}
