/* output.c -- the output field table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "output.h"

const out_field OUT_FIELDS[OUT_N_FIELDS] = {
    [OUT_COVERAGE]   = { "coverage",       ACCUM_COVERAGE,  false, OUT_ADD       },
    [OUT_REACTIVITY] = { "reactivity",     ACCUM_MUTATIONS, false, OUT_SUBTRACT  },
    [OUT_ERROR]      = { "error",          ACCUM_MUTATIONS, false, OUT_PROPAGATE },
    [OUT_LENGTHS]    = { "reads/lengths",  ACCUM_LENGTHS,   true,  OUT_ADD       },
    [OUT_READS]      = { "reads/counted",  ACCUM_READS,     true,  OUT_ADD       },
    [OUT_REJECTED]   = { "reads/rejected", ACCUM_FILTERED,  true,  OUT_ADD       },
};

size_t out_extent(out_field_id id, size_t len, size_t cap)
{
    return accum_extent(OUT_FIELDS[id].shape, len, cap);
}

size_t out_widest(size_t cap)
{
    size_t widest = 0;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        size_t width = out_extent(id, cap, cap);

        widest = width > widest ? width : widest;
    }

    return widest;
}

int out_rank(out_field_id id)
{
    return ACCUM_FIELDS[OUT_FIELDS[id].shape].kind == ACCUM_SCALAR
         ? OUT_RANK_SCALAR : OUT_RANK_VECTOR;
}
