/* h5layout.c -- datatypes, shapes, dataspaces and property lists.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "h5layout.h"

#include <math.h>
#include <stdbool.h>

/* Chunks are sized in bytes and not rows, so that a file of few long references
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

static hsize_t rows_per_chunk(size_t row_bytes, int32_t n_refs)
{
    hsize_t rows = row_bytes ? TARGET_CHUNK_BYTES / row_bytes : (hsize_t)n_refs;

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
    size_t extents[OUT_RANK_MAX];
    int    rank = out_dims(id, n_refs, cap, extents);

    for (int i = 0; i < rank; i++) {
        dims[i] = chunk[i] = (hsize_t)extents[i];
    }

    /* A chunk spans whole rows and as many of them as fit, so only the reference
     * dimension is cut down. A field with no reference dimension is one value, and is
     * not chunked at all. */
    if (OUT_FIELDS[id].per_ref) {
        chunk[0] = rows_per_chunk(out_values(id, cap, cap) * out_stored_bytes(id), n_refs);
    }
}

/* ------------------------------------------------------------------------ */
/* Storage                                                                   */
/* ------------------------------------------------------------------------ */

hid_t h5layout_type(out_field_id id)
{
    switch (OUT_FIELDS[id].stored) {
        case OUT_F32:      return H5T_IEEE_F32LE;
        case OUT_U64:      return H5T_STD_U64LE;
        case OUT_N_STORED: break;
    }

    return H5I_INVALID_HID;
}

/* Returns the type a field's values are handed over in, which is the one it is stored
 * as. A value read from a file and written to another passes through nothing wider on
 * the way. */
hid_t h5layout_memory_type(out_field_id id)
{
    switch (OUT_FIELDS[id].stored) {
        case OUT_F32:      return H5T_NATIVE_FLOAT;
        case OUT_U64:      return H5T_NATIVE_UINT64;
        case OUT_N_STORED: break;
    }

    return H5I_INVALID_HID;
}

/* Gives a field's absence marker, in the type the field is stored as.
 *
 * An unsigned has no NaN, so that one pairing has no fill and a field declaring it is
 * refused, and not filled with something that means otherwise. A storage type added
 * without fills of its own is refused the same way. */
const void *h5layout_fill(out_field_id id)
{
    static const uint64_t u64_zero = 0;
    static const float    f32_zero = 0.0f;
    static const float    f32_nan  = (float)NAN;

    static const void *const fill[OUT_N_STORED][OUT_N_ABSENT] = {
        [OUT_F32] = { [OUT_ZERO] = &f32_zero, [OUT_NAN] = &f32_nan },
        [OUT_U64] = { [OUT_ZERO] = &u64_zero },
    };

    return fill[OUT_FIELDS[id].stored][OUT_FIELDS[id].absent];
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
 * it: the arrays are filled in full and the rank bounds how much is read. */
int h5layout_select_span(hid_t filespace, hid_t memspace, out_field_id id,
                         int32_t tid, size_t from, size_t n)
{
    hsize_t start[OUT_RANK_MAX] = { (hsize_t)tid, (hsize_t)from };
    hsize_t count[OUT_RANK_MAX] = { 1, (hsize_t)n };
    hsize_t offset              = 0;
    hsize_t extent              = out_rank(id) > 1 ? (hsize_t)n : 1;

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
    const void *fill = h5layout_fill(id);
    hid_t       dcpl = h5layout_untimed_plist(H5P_DATASET_CREATE);

    if (dcpl < 0) {
        return H5I_INVALID_HID;
    }

    /* A dataset of rank zero holds one value, which is written the moment it is
     * created. There is nothing to chunk, nothing to filter, and no position that
     * could go unwritten and need a fill. */
    if (rank == 0) {
        return dcpl;
    }

    if (!fill ||
        H5Pset_chunk(dcpl, rank, chunk) < 0 ||
        H5Pset_fill_value(dcpl, h5layout_type(id), fill) < 0 ||
        H5Pset_fill_time(dcpl, H5D_FILL_TIME_ALLOC) < 0 ||
        H5Pset_shuffle(dcpl) < 0 ||
        H5Pset_deflate(dcpl, DEFLATE_LEVEL) < 0) {
        H5Pclose(dcpl);
        return H5I_INVALID_HID;
    }

    return dcpl;
}

hid_t h5layout_access_plist(out_field_id id, const hsize_t *chunk, int rank)
{
    hid_t  dapl  = H5Pcreate(H5P_DATASET_ACCESS);
    size_t bytes = out_stored_bytes(id);

    if (dapl < 0) {
        return H5I_INVALID_HID;
    }

    /* Nothing to cache where nothing is chunked. */
    if (rank == 0) {
        return dapl;
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
