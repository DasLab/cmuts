/* h5reader.c -- HDF5 input, one row per reference.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "h5reader.h"

#include <stdio.h>
#include <stdlib.h>

#include "error.h"
#include "h5layout.h"

/* The field whose shape gives the rest. Coverage is one value per base, so the width
 * of its row is the capacity every other field's width derives from. */
#define SHAPE_FIELD OUT_COVERAGE

struct h5reader {
    hid_t   file;
    hid_t   dataset[OUT_N_FIELDS];
    /* A dataspace apiece, kept for the life of the reader as the writer keeps its
     * own: only the selection differs between one row and the next. */
    hid_t   filespace[OUT_N_FIELDS];
    hid_t   memspace;   /* one row of the widest field, selected down to size */
    int32_t n_refs;
    size_t  ref_cap;
    char    error[CM_ERROR_MAX];
};

static int fail(h5reader *r, const char *what)
{
    snprintf(r->error, sizeof r->error, "%s", what);
    return -1;
}

static int fail_field(h5reader *r, out_field_id id, const char *what)
{
    snprintf(r->error, sizeof r->error, "%s: %s", OUT_FIELDS[id].name, what);
    return -1;
}

/* ------------------------------------------------------------------------ */
/* Shape                                                                     */
/* ------------------------------------------------------------------------ */

/* A dataset's dimensions, failing where its rank is not the one expected. */
static int dataset_dims(hid_t dataset, int rank, hsize_t *dims)
{
    hid_t space = H5Dget_space(dataset);
    int   found;

    if (space < 0) {
        return -1;
    }

    found = H5Sget_simple_extent_ndims(space);

    if (found != rank || H5Sget_simple_extent_dims(space, dims, NULL) < 0) {
        H5Sclose(space);
        return -1;
    }

    H5Sclose(space);
    return 0;
}

/* Reads the number of references and the capacity from SHAPE_FIELD, whose width is
 * the capacity itself. Every other field is then checked against them. */
static int probe_shape(h5reader *r)
{
    hid_t   dataset = H5Dopen2(r->file, OUT_FIELDS[SHAPE_FIELD].name, H5P_DEFAULT);
    hsize_t dims[OUT_RANK_MAX];
    int     status;

    if (dataset < 0) {
        return fail_field(r, SHAPE_FIELD, "not present; this is not a cmuts output");
    }

    status = dataset_dims(dataset, out_rank(SHAPE_FIELD), dims);
    H5Dclose(dataset);

    if (status < 0) {
        return fail_field(r, SHAPE_FIELD, "is not a row per reference");
    }

    if (dims[0] == 0 || dims[0] > INT32_MAX || dims[1] == 0) {
        return fail(r, "the file holds no references");
    }

    r->n_refs  = (int32_t)dims[0];
    r->ref_cap = (size_t)dims[1];

    return 0;
}

static int check_shape(h5reader *r, out_field_id id, const hsize_t *expected)
{
    int     rank = out_rank(id);
    hsize_t dims[OUT_RANK_MAX];

    if (dataset_dims(r->dataset[id], rank, dims) < 0) {
        return fail_field(r, id, "has an unexpected number of dimensions");
    }

    for (int i = 0; i < rank; i++) {
        if (dims[i] != expected[i]) {
            return fail_field(r, id, "does not agree with the shape of the file");
        }
    }

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Lifetime                                                                  */
/* ------------------------------------------------------------------------ */

static int open_field(h5reader *r, out_field_id id)
{
    hsize_t dims[OUT_RANK_MAX]  = { 0, 0 };
    hsize_t chunk[OUT_RANK_MAX] = { 0, 0 };
    hid_t   dapl;

    h5layout_shape(id, r->n_refs, r->ref_cap, dims, chunk);

    dapl = h5layout_access_plist(id, chunk, out_rank(id));
    if (dapl < 0) {
        return fail(r, "unable to prepare a dataset for reading");
    }

    r->dataset[id] = H5Dopen2(r->file, OUT_FIELDS[id].name, dapl);
    H5Pclose(dapl);

    if (r->dataset[id] < 0) {
        return fail_field(r, id, "not present");
    }

    if (check_shape(r, id, dims) < 0) {
        return -1;
    }

    r->filespace[id] = H5Dget_space(r->dataset[id]);

    return r->filespace[id] < 0 ? fail_field(r, id, "cannot be described") : 0;
}

static int open_fields(h5reader *r)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (open_field(r, id) < 0) {
            return -1;
        }
    }

    return 0;
}

/* Allocates a reader holding nothing yet, with every handle marked absent.
 *
 * The steps that build the rest may each fail and leave those after them undone, and the
 * reader is closed whatever happened, so it must be safe to close from here onwards: it
 * closes exactly what it opened. Zero, which calloc leaves behind, is a handle HDF5 would
 * accept, hence the marking. */
static h5reader *reader_alloc(void)
{
    h5reader *r = calloc(1, sizeof *r);

    if (!r) {
        return NULL;
    }

    /* Report failures through h5reader_error rather than HDF5's own stack trace on
     * stderr. */
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        r->dataset[id]   = H5I_INVALID_HID;
        r->filespace[id] = H5I_INVALID_HID;
    }

    r->file     = H5I_INVALID_HID;
    r->memspace = H5I_INVALID_HID;

    return r;
}

static int open_file(h5reader *r, const char *path)
{
    r->file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);

    return r->file < 0
         ? fail(r, "unable to open the file, which may not be HDF5 at all") : 0;
}

/* Prepares the row every read is selected into, which the shape must be known to size. */
static int build_memspace(h5reader *r)
{
    r->memspace = h5layout_row_space(r->ref_cap);

    return r->memspace < 0 ? fail(r, "unable to prepare the file for reading") : 0;
}

h5reader *h5reader_open(const char *path)
{
    h5reader *r = reader_alloc();

    if (!r) {
        return NULL;
    }

    if (open_file(r, path) == 0 &&
        probe_shape(r) == 0 &&
        build_memspace(r) == 0) {
        open_fields(r);
    }

    return r;
}

void h5reader_close(h5reader *r)
{
    if (!r) {
        return;
    }

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (r->filespace[id] >= 0) {
            H5Sclose(r->filespace[id]);
        }
        if (r->dataset[id] >= 0) {
            H5Dclose(r->dataset[id]);
        }
    }

    if (r->memspace >= 0) {
        H5Sclose(r->memspace);
    }

    if (r->file >= 0) {
        H5Fclose(r->file);
    }

    free(r);
}

const char *h5reader_error(const h5reader *r)
{
    return r->error[0] ? r->error : NULL;
}

int32_t h5reader_refs(const h5reader *r)
{
    return r->n_refs;
}

size_t h5reader_capacity(const h5reader *r)
{
    return r->ref_cap;
}

/* ------------------------------------------------------------------------ */
/* Rows                                                                      */
/* ------------------------------------------------------------------------ */

int h5reader_field(h5reader *r, out_field_id id, int32_t tid, double *values)
{
    size_t width = out_values(id, r->ref_cap, r->ref_cap);
    herr_t status;

    if (tid < 0 || tid >= r->n_refs) {
        return fail(r, "reference index outside the file");
    }

    if (h5layout_select_span(r->filespace[id], r->memspace, id, tid, 0, width) < 0) {
        return fail(r, "unable to select an input row");
    }

    status = H5Dread(r->dataset[id], H5T_NATIVE_DOUBLE, r->memspace,
                     r->filespace[id], H5P_DEFAULT, values);

    return status < 0 ? fail(r, "unable to read an input row") : 0;
}

/* ------------------------------------------------------------------------ */
/* Totals                                                                    */
/* ------------------------------------------------------------------------ */

int h5reader_total(h5reader *r, out_field_id id, size_t *value)
{
    uint64_t stored = 0;
    herr_t   status;

    if (OUT_FIELDS[id].per_ref) {
        return fail(r, "a field with a row per reference has no run total");
    }

    status = H5Dread(r->dataset[id], H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, &stored);

    if (status < 0) {
        return fail_field(r, id, "unable to read it");
    }

    *value = (size_t)stored;
    return 0;
}
