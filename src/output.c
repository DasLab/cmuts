/* output.c -- the output field table and the HDF5 property lists it defines.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "output.h"

#include <math.h>

/* Chunks are sized by bytes rather than rows so that a file of few long
 * references and one of many short references both land near this figure. */
#define TARGET_CHUNK_BYTES (1u << 20)

/* The slots are a hash table over chunk indices, which is why their number is
 * prime. */
#define CACHED_CHUNKS 4
#define CACHE_SLOTS   521

/* HDF5's own default, restated because naming either of the other two means
 * passing all three. */
#define CACHE_PREEMPTION 0.75

#define DEFLATE_LEVEL 3

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

/* ------------------------------------------------------------------------ */
/* Shape                                                                     */
/* ------------------------------------------------------------------------ */

static hsize_t rows_per_chunk(size_t row_values, int32_t n_refs)
{
    size_t  row_bytes = row_values * sizeof(float);
    hsize_t rows      = row_bytes ? TARGET_CHUNK_BYTES / row_bytes : (hsize_t)n_refs;

    if (rows < 1)
        rows = 1;
    if (rows > (hsize_t)n_refs)
        rows = (hsize_t)n_refs;

    return rows;
}

void out_shape(out_field_id id, int32_t n_refs, size_t cap,
               hsize_t *dims, hsize_t *chunk)
{
    size_t width = out_extent(id, cap, cap);

    dims[0]  = (hsize_t)n_refs;
    chunk[0] = rows_per_chunk(width, n_refs);

    if (out_rank(id) == OUT_RANK_VECTOR) {
        dims[1]  = width;
        chunk[1] = width;
    }
}

/* ------------------------------------------------------------------------ */
/* Storage                                                                   */
/* ------------------------------------------------------------------------ */

hid_t out_type(out_field_id id)
{
    return OUT_FIELDS[id].counted ? H5T_STD_U64LE : H5T_IEEE_F32LE;
}

/* Zero for a count and for what is measured, NaN for a rate that was not.
 *
 * A position no read reached was reached by no read, which is what a count of
 * zero says. A rate is not a count: a reference no read named has no rate, and
 * filling one with zero would say its every position was measured and found
 * unmodified, which is the most confident thing the output can say and it would
 * be saying it about nothing at all. A count has no NaN to be had, being an
 * unsigned, and needs none: nothing it is written for has padding. */
const void *out_fill(out_field_id id)
{
    static const uint64_t none = 0;
    static const float    zero = 0.0f;
    static const float    nan  = (float)NAN;

    if (OUT_FIELDS[id].counted)
        return &none;

    return id == OUT_REACTIVITY || id == OUT_ERROR ? (const void *)&nan
                                                   : (const void *)&zero;
}

hid_t out_untimed_plist(hid_t class_id)
{
    hid_t plist = H5Pcreate(class_id);

    if (plist < 0)
        return H5I_INVALID_HID;

    if (H5Pset_obj_track_times(plist, false) < 0) {
        H5Pclose(plist);
        return H5I_INVALID_HID;
    }

    return plist;
}

hid_t out_layout_plist(out_field_id id, const hsize_t *chunk, int rank)
{
    hid_t dcpl = out_untimed_plist(H5P_DATASET_CREATE);

    if (dcpl < 0)
        return H5I_INVALID_HID;

    if (H5Pset_chunk(dcpl, rank, chunk) < 0 ||
        H5Pset_fill_value(dcpl, out_type(id), out_fill(id)) < 0 ||
        H5Pset_fill_time(dcpl, H5D_FILL_TIME_ALLOC) < 0 ||
        H5Pset_shuffle(dcpl) < 0 ||
        H5Pset_deflate(dcpl, DEFLATE_LEVEL) < 0) {
        H5Pclose(dcpl);
        return H5I_INVALID_HID;
    }

    return dcpl;
}

hid_t out_access_plist(const hsize_t *chunk, int rank)
{
    hid_t  dapl  = H5Pcreate(H5P_DATASET_ACCESS);
    size_t bytes = sizeof(float);

    if (dapl < 0)
        return H5I_INVALID_HID;

    for (int i = 0; i < rank; i++)
        bytes *= (size_t)chunk[i];

    if (H5Pset_chunk_cache(dapl, CACHE_SLOTS, CACHED_CHUNKS * bytes,
                           CACHE_PREEMPTION) < 0) {
        H5Pclose(dapl);
        return H5I_INVALID_HID;
    }

    return dapl;
}
