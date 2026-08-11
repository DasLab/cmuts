/* h5layout.c -- datatypes, shapes, dataspaces and property lists.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "h5layout.h"

#include <math.h>

/* Chunks are sized in bytes rather than rows, so that a file of few long references
 * and one of many short references both land near this figure. */
#define TARGET_CHUNK_BYTES (1u << 20)

/* The slots are a hash table over chunk indices, hence a prime. */
#define CACHED_CHUNKS 4
#define CACHE_SLOTS   521

/* HDF5's own default, restated because setting either of the other two means
 * passing all three. */
#define CACHE_PREEMPTION 0.75

#define DEFLATE_LEVEL 3

/* ------------------------------------------------------------------------ */
/* Shape                                                                     */
/* ------------------------------------------------------------------------ */

static hsize_t rows_per_chunk(size_t row_values, int32_t n_refs)
{
    size_t  row_bytes = row_values * sizeof(float);
    hsize_t rows      = row_bytes ? TARGET_CHUNK_BYTES / row_bytes : (hsize_t)n_refs;

    if (rows < 1) {
        rows = 1;
    }
    if (rows > (hsize_t)n_refs) {
        rows = (hsize_t)n_refs;
    }

    return rows;
}

void h5layout_shape(out_field_id id, int32_t n_refs, size_t cap,
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

hid_t h5layout_type(out_field_id id)
{
    return OUT_FIELDS[id].counted ? H5T_STD_U64LE : H5T_IEEE_F32LE;
}

/* A field's OUT_ZERO or OUT_NAN, in the type the field is stored as. An unsigned has
 * no NaN, so OUT_NAN is not available to a counted field; the table gives each of
 * them OUT_ZERO, no counted field having padding. */
const void *h5layout_fill(out_field_id id)
{
    static const uint64_t none = 0;
    static const float    zero = 0.0f;
    static const float    nan  = (float)NAN;

    if (OUT_FIELDS[id].counted) {
        return &none;
    }

    return OUT_FIELDS[id].absent == OUT_NAN ? (const void *)&nan
                                            : (const void *)&zero;
}

/* ------------------------------------------------------------------------ */
/* Dataspaces                                                                */
/* ------------------------------------------------------------------------ */

hid_t h5layout_row_space(size_t cap)
{
    hsize_t widest = (hsize_t)out_widest(cap);

    return H5Screate_simple(1, &widest, NULL);
}

/* Selects n values of one reference's row, in the file and in memory together. A
 * scalar field has one value per reference, so the column plays no part in selecting
 * it: the arrays are filled in full and the rank decides how much is read. */
int h5layout_select_span(hid_t filespace, hid_t memspace, out_field_id id,
                         int32_t tid, size_t from, size_t n)
{
    hsize_t start[OUT_RANK_MAX] = { (hsize_t)tid, (hsize_t)from };
    hsize_t count[OUT_RANK_MAX] = { 1, (hsize_t)n };
    hsize_t offset              = 0;
    hsize_t extent              = out_rank(id) == OUT_RANK_VECTOR ? (hsize_t)n : 1;

    if (H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
        return -1;
    }

    return H5Sselect_hyperslab(memspace, H5S_SELECT_SET, &offset, NULL,
                               &extent, NULL);
}

/* ------------------------------------------------------------------------ */
/* Property lists                                                            */
/* ------------------------------------------------------------------------ */

hid_t h5layout_untimed_plist(hid_t class_id)
{
    hid_t plist = H5Pcreate(class_id);

    if (plist < 0) {
        return H5I_INVALID_HID;
    }

    if (H5Pset_obj_track_times(plist, false) < 0) {
        H5Pclose(plist);
        return H5I_INVALID_HID;
    }

    return plist;
}

hid_t h5layout_creation_plist(out_field_id id, const hsize_t *chunk, int rank)
{
    hid_t dcpl = h5layout_untimed_plist(H5P_DATASET_CREATE);

    if (dcpl < 0) {
        return H5I_INVALID_HID;
    }

    if (H5Pset_chunk(dcpl, rank, chunk) < 0 ||
        H5Pset_fill_value(dcpl, h5layout_type(id), h5layout_fill(id)) < 0 ||
        H5Pset_fill_time(dcpl, H5D_FILL_TIME_ALLOC) < 0 ||
        H5Pset_shuffle(dcpl) < 0 ||
        H5Pset_deflate(dcpl, DEFLATE_LEVEL) < 0) {
        H5Pclose(dcpl);
        return H5I_INVALID_HID;
    }

    return dcpl;
}

hid_t h5layout_access_plist(const hsize_t *chunk, int rank)
{
    hid_t  dapl  = H5Pcreate(H5P_DATASET_ACCESS);
    size_t bytes = sizeof(float);

    if (dapl < 0) {
        return H5I_INVALID_HID;
    }

    for (int i = 0; i < rank; i++) {
        bytes *= (size_t)chunk[i];
    }

    if (H5Pset_chunk_cache(dapl, CACHE_SLOTS, CACHED_CHUNKS * bytes,
                           CACHE_PREEMPTION) < 0) {
        H5Pclose(dapl);
        return H5I_INVALID_HID;
    }

    return dapl;
}
