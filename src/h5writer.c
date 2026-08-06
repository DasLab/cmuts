/* h5writer.c -- HDF5 output, one row per reference.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "h5writer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <hdf5.h>

#define H5WRITER_ERROR_MAX 512

#define DATASET_REFERENCE "reference"

/* Chunks are sized by bytes rather than rows so that a file of few long
 * references and one of many short references both land near this figure. */
#define TARGET_CHUNK_BYTES (1u << 20)

/* Rows are completed in roughly ascending order, so holding a handful of
 * chunks in the cache keeps a single-row write from forcing a read, modify and
 * rewrite of the chunk around it. */
#define CACHED_CHUNKS 4
#define CACHE_SLOTS   521

#define DEFLATE_LEVEL 4

struct h5writer {
    hid_t   file;
    hid_t   dataset[ACCUM_N_FIELDS];
    int32_t n_refs;
    size_t  ref_cap;
    char    error[H5WRITER_ERROR_MAX];
};

static int fail(h5writer *w, const char *what)
{
    snprintf(w->error, sizeof w->error, "%s", what);
    return -1;
}

/* ------------------------------------------------------------------------ */
/* Dataset construction                                                      */
/* ------------------------------------------------------------------------ */

static hsize_t rows_per_chunk(const h5writer *w)
{
    size_t  row_bytes = w->ref_cap * sizeof(float);
    hsize_t rows      = row_bytes ? TARGET_CHUNK_BYTES / row_bytes : (hsize_t)w->n_refs;

    if (rows < 1)
        rows = 1;
    if (rows > (hsize_t)w->n_refs)
        rows = (hsize_t)w->n_refs;

    return rows;
}

static int field_rank(accum_field_id id)
{
    return ACCUM_FIELDS[id].kind == ACCUM_PER_BASE ? 2 : 1;
}

static void field_shape(const h5writer *w, accum_field_id id, hsize_t *dims, hsize_t *chunk)
{
    dims[0]  = (hsize_t)w->n_refs;
    chunk[0] = rows_per_chunk(w);

    if (field_rank(id) == 2) {
        dims[1]  = w->ref_cap;
        chunk[1] = w->ref_cap;
    }
}

static hid_t make_layout(const hsize_t *chunk, int rank)
{
    hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
    float fill = (float)NAN;

    if (dcpl < 0)
        return H5I_INVALID_HID;

    if (H5Pset_chunk(dcpl, rank, chunk) < 0 ||
        H5Pset_fill_value(dcpl, H5T_NATIVE_FLOAT, &fill) < 0 ||
        H5Pset_fill_time(dcpl, H5D_FILL_TIME_ALLOC) < 0 ||
        H5Pset_shuffle(dcpl) < 0 ||
        H5Pset_deflate(dcpl, DEFLATE_LEVEL) < 0) {
        H5Pclose(dcpl);
        return H5I_INVALID_HID;
    }

    return dcpl;
}

static hid_t make_access(const hsize_t *chunk, int rank)
{
    hid_t  dapl  = H5Pcreate(H5P_DATASET_ACCESS);
    size_t bytes = sizeof(float);

    if (dapl < 0)
        return H5I_INVALID_HID;

    for (int i = 0; i < rank; i++)
        bytes *= (size_t)chunk[i];

    if (H5Pset_chunk_cache(dapl, CACHE_SLOTS, CACHED_CHUNKS * bytes, 0.75) < 0) {
        H5Pclose(dapl);
        return H5I_INVALID_HID;
    }

    return dapl;
}

static hid_t create_field(h5writer *w, accum_field_id id)
{
    hsize_t dims[2]  = { 0, 0 };
    hsize_t chunk[2] = { 0, 0 };
    int     rank     = field_rank(id);
    hid_t   space, dcpl, dapl, dataset;

    field_shape(w, id, dims, chunk);

    space = H5Screate_simple(rank, dims, NULL);
    dcpl  = make_layout(chunk, rank);
    dapl  = make_access(chunk, rank);

    if (space < 0 || dcpl < 0 || dapl < 0) {
        dataset = H5I_INVALID_HID;
    } else {
        dataset = H5Dcreate2(w->file, ACCUM_FIELDS[id].name, H5T_IEEE_F32LE,
                             space, H5P_DEFAULT, dcpl, dapl);
    }

    H5Pclose(dapl);
    H5Pclose(dcpl);
    H5Sclose(space);

    return dataset;
}

/* ------------------------------------------------------------------------ */
/* Lifetime                                                                  */
/* ------------------------------------------------------------------------ */

h5writer *h5writer_create(const char *path, int32_t n_refs, size_t ref_cap,
                          bool overwrite)
{
    h5writer *w = calloc(1, sizeof *w);
    if (!w)
        return NULL;

    /* Report failures through h5writer_error rather than HDF5's own stack
     * trace on stderr. */
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    w->n_refs  = n_refs;
    w->ref_cap = ref_cap;

    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++)
        w->dataset[id] = H5I_INVALID_HID;

    /* Exclusive unless replacing was asked for, so that the decision cannot
     * be undone by the file appearing between the check and the create. */
    w->file = H5Fcreate(path, overwrite ? H5F_ACC_TRUNC : H5F_ACC_EXCL,
                        H5P_DEFAULT, H5P_DEFAULT);
    if (w->file < 0) {
        fail(w, "unable to create the output file");
        return w;
    }

    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++) {
        w->dataset[id] = create_field(w, id);
        if (w->dataset[id] < 0) {
            fail(w, "unable to create a dataset");
            return w;
        }
    }

    return w;
}

void h5writer_close(h5writer *w)
{
    if (!w)
        return;

    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++)
        if (w->dataset[id] >= 0)
            H5Dclose(w->dataset[id]);

    if (w->file >= 0)
        H5Fclose(w->file);

    free(w);
}

const char *h5writer_error(const h5writer *w)
{
    return w->error[0] ? w->error : NULL;
}

/* ------------------------------------------------------------------------ */
/* Rows                                                                      */
/* ------------------------------------------------------------------------ */

/* Selects the part of a dataset belonging to one reference: the whole row for
 * a scalar field, its first len values for a per-base field. */
static int select_row(hid_t dataset, int32_t tid, size_t len, int rank,
                      hid_t *filespace, hid_t *memspace)
{
    hsize_t start[2] = { (hsize_t)tid, 0 };
    hsize_t count[2] = { 1, (hsize_t)len };
    hsize_t extent   = rank == 2 ? (hsize_t)len : 1;

    *filespace = H5Dget_space(dataset);
    if (*filespace < 0)
        return -1;

    if (H5Sselect_hyperslab(*filespace, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
        H5Sclose(*filespace);
        return -1;
    }

    *memspace = H5Screate_simple(1, &extent, NULL);
    if (*memspace < 0) {
        H5Sclose(*filespace);
        return -1;
    }

    return 0;
}

static int write_field(h5writer *w, accum_field_id id, int32_t tid, size_t len,
                       const accum *acc)
{
    int   rank = field_rank(id);
    hid_t filespace, memspace;
    herr_t status;

    if (select_row(w->dataset[id], tid, len, rank, &filespace, &memspace) < 0)
        return fail(w, "unable to select an output row");

    status = H5Dwrite(w->dataset[id], H5T_NATIVE_DOUBLE, memspace, filespace,
                      H5P_DEFAULT, accum_const_data(acc, id));

    H5Sclose(memspace);
    H5Sclose(filespace);

    return status < 0 ? fail(w, "unable to write an output row") : 0;
}

int h5writer_row(h5writer *w, int32_t tid, size_t len, const accum *acc)
{
    if (tid < 0 || tid >= w->n_refs)
        return fail(w, "reference index outside the output");

    if (len == 0)
        return 0;

    for (accum_field_id id = 0; id < ACCUM_N_FIELDS; id++)
        if (write_field(w, id, tid, len, acc) < 0)
            return -1;

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Names and totals                                                          */
/* ------------------------------------------------------------------------ */

static hid_t make_string_type(void)
{
    hid_t type = H5Tcopy(H5T_C_S1);

    if (type < 0)
        return H5I_INVALID_HID;

    if (H5Tset_size(type, H5T_VARIABLE) < 0 || H5Tset_cset(type, H5T_CSET_UTF8) < 0) {
        H5Tclose(type);
        return H5I_INVALID_HID;
    }

    return type;
}

int h5writer_names(h5writer *w, const char *const *names, int32_t n_refs)
{
    hsize_t dims    = (hsize_t)n_refs;
    hid_t   type    = make_string_type();
    hid_t   space   = H5Screate_simple(1, &dims, NULL);
    hid_t   dataset = H5I_INVALID_HID;
    herr_t  status  = -1;

    if (type >= 0 && space >= 0)
        dataset = H5Dcreate2(w->file, DATASET_REFERENCE, type, space,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (dataset >= 0)
        status = H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, names);

    if (dataset >= 0)
        H5Dclose(dataset);
    if (space >= 0)
        H5Sclose(space);
    if (type >= 0)
        H5Tclose(type);

    return status < 0 ? fail(w, "unable to write the reference names") : 0;
}

int h5writer_count(h5writer *w, const char *name, size_t value)
{
    uint64_t stored  = value;
    hid_t    space   = H5Screate(H5S_SCALAR);
    hid_t    attr    = H5I_INVALID_HID;
    herr_t   status  = -1;

    if (space >= 0)
        attr = H5Acreate2(w->file, name, H5T_NATIVE_UINT64, space,
                          H5P_DEFAULT, H5P_DEFAULT);

    if (attr >= 0)
        status = H5Awrite(attr, H5T_NATIVE_UINT64, &stored);

    if (attr >= 0)
        H5Aclose(attr);
    if (space >= 0)
        H5Sclose(space);

    return status < 0 ? fail(w, "unable to write a run total") : 0;
}
