/* output.c -- the output field table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "output.h"

const out_field OUT_FIELDS[OUT_N_FIELDS] = {
    [OUT_COVERAGE]   = { "coverage",       SHAPE_PER_BASE,   OUT_F32, OUT_ZERO, OUT_ADD       },
    [OUT_REACTIVITY] = { "reactivity",     SHAPE_PER_BASE,   OUT_F32, OUT_NAN,  OUT_SUBTRACT  },
    [OUT_ERROR]      = { "error",          SHAPE_PER_BASE,   OUT_F32, OUT_NAN,  OUT_PROPAGATE },
    [OUT_LENGTHS]    = { "reads/lengths",  SHAPE_PER_LENGTH, OUT_U64, OUT_ZERO, OUT_ADD       },
    [OUT_READS]      = { "reads/counted",  SHAPE_SCALAR,     OUT_U64, OUT_ZERO, OUT_ADD       },
    [OUT_REJECTED]   = { "reads/rejected", SHAPE_SCALAR,     OUT_U64, OUT_ZERO, OUT_ADD       },
};

size_t out_extent(out_field_id id, size_t len, size_t cap)
{
    return shape_extent(OUT_FIELDS[id].shape, len, cap);
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
    return shape_rank(OUT_FIELDS[id].shape);
}
