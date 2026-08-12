/* output.c -- the output field table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "output.h"

const out_field OUT_FIELDS[OUT_N_FIELDS] = {
    [OUT_COVERAGE]   = { "coverage",       shape_per_base,   true,  OUT_F32, OUT_ZERO },
    [OUT_REACTIVITY] = { "reactivity",     shape_per_base,   true,  OUT_F32, OUT_NAN  },
    [OUT_ERROR]      = { "error",          shape_per_base,   true,  OUT_F32, OUT_NAN  },
    [OUT_LENGTHS]    = { "reads/lengths",  shape_per_length, true,  OUT_U64, OUT_ZERO },
    [OUT_READS]      = { "reads/counted",  shape_none,       true,  OUT_U64, OUT_ZERO },
    [OUT_REJECTED]   = { "reads/rejected", shape_none,       true,  OUT_U64, OUT_ZERO },
    [OUT_UNMAPPED]   = { "reads/unmapped", shape_none,       false, OUT_U64, OUT_ZERO },
};

size_t out_values(out_field_id id, size_t len, size_t cap)
{
    return shape_values(OUT_FIELDS[id].row, len, cap);
}

size_t out_widest(size_t cap)
{
    size_t widest = 0;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        size_t width = out_values(id, cap, cap);

        widest = width > widest ? width : widest;
    }

    return widest;
}

size_t out_stored_bytes(out_field_id id)
{
    switch (OUT_FIELDS[id].stored) {
        case OUT_F32:      return sizeof(float);
        case OUT_U64:      return sizeof(uint64_t);
        case OUT_N_STORED: break;
    }

    return 0;
}

size_t out_widest_bytes(void)
{
    size_t widest = 0;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        size_t bytes = out_stored_bytes(id);

        widest = bytes > widest ? bytes : widest;
    }

    return widest;
}

int out_dims(out_field_id id, int32_t n_refs, size_t cap, size_t *dims)
{
    shape_extents row     = OUT_FIELDS[id].row(cap, cap);
    int           extents = shape_rank(row);
    int           rank    = 0;

    if (OUT_FIELDS[id].per_ref) {
        dims[rank++] = (size_t)n_refs;
    }

    for (int i = 0; i < extents; i++) {
        dims[rank++] = row.dim[i];
    }

    return rank;
}

/* Gives the rank, discarding the dimensions an empty run produces. */
int out_rank(out_field_id id)
{
    size_t dims[OUT_RANK_MAX];

    return out_dims(id, 0, 0, dims);
}
