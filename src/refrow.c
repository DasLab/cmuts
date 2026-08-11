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

    if (!r)
        return NULL;

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
    if (!r)
        return;

    free(r->row);
    free(r);
}

/* One field's values. The reactivity and its error are derived rather than accumulated,
 * so they are computed into the scratch row; every other field is read in place. */
static const double *values(refrow *r, out_field_id id, const accum *acc,
                            size_t len)
{
    if (id == OUT_REACTIVITY)
        rate_reactivity(&r->rates, acc, len, r->row);
    else if (id == OUT_ERROR)
        rate_error(&r->rates, acc, len, r->row);
    else
        return accum_const_data(acc, OUT_FIELDS[id].shape);

    return r->row;
}

int refrow_write(refrow *r, int32_t tid, size_t len, const accum *acc)
{
    if (len == 0)
        return 0;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++)
        if (h5writer_field(r->out, id, tid, len, values(r, id, acc, len)) < 0)
            return -1;

    return 0;
}
