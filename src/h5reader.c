/* h5reader.c -- HDF5 input, one row per reference.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "h5reader.h"

#include <stdio.h>
#include <stdlib.h>

#include "error.h"
#include "output.h"

/* The field whose shape gives the rest. Coverage is one value per base, so the
 * width of its row is the capacity every other field's width derives from. */
#define SHAPE_FIELD OUT_COVERAGE

struct h5reader {
    hid_t   file;
    hid_t   dataset[OUT_N_FIELDS];
    hid_t   reads;      /* the group the per-run counts are gathered in */
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

/* A dataset's own dimensions, which are checked against what the layout says
 * they should be rather than trusted. */
static int dataset_dims(hid_t dataset, int rank, hsize_t *dims)
{
    hid_t space = H5Dget_space(dataset);
    int   found;

    if (space < 0)
        return -1;

    found = H5Sget_simple_extent_ndims(space);

    if (found != rank || H5Sget_simple_extent_dims(space, dims, NULL) < 0) {
        H5Sclose(space);
        return -1;
    }

    H5Sclose(space);
    return 0;
}

/* The number of references and the capacity, taken from the one field whose
 * width is the capacity itself. Every other field is then checked against them.
 */
static int probe_shape(h5reader *r)
{
    hid_t   dataset = H5Dopen2(r->file, OUT_FIELDS[SHAPE_FIELD].name, H5P_DEFAULT);
    hsize_t dims[OUT_RANK_MAX];
    int     status;

    if (dataset < 0)
        return fail_field(r, SHAPE_FIELD, "not present; this is not a cmuts output");

    status = dataset_dims(dataset, OUT_RANK_VECTOR, dims);
    H5Dclose(dataset);

    if (status < 0)
        return fail_field(r, SHAPE_FIELD, "is not a row per reference");

    if (dims[0] == 0 || dims[0] > INT32_MAX || dims[1] == 0)
        return fail(r, "the file holds no references");

    r->n_refs  = (int32_t)dims[0];
    r->ref_cap = (size_t)dims[1];

    return 0;
}

static int check_shape(h5reader *r, out_field_id id, const hsize_t *expected)
{
    int     rank = out_rank(id);
    hsize_t dims[OUT_RANK_MAX];

    if (dataset_dims(r->dataset[id], rank, dims) < 0)
        return fail_field(r, id, "has an unexpected number of dimensions");

    for (int i = 0; i < rank; i++)
        if (dims[i] != expected[i])
            return fail_field(r, id, "does not agree with the shape of the file");

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

    out_shape(id, r->n_refs, r->ref_cap, dims, chunk);

    dapl = out_access_plist(chunk, out_rank(id));
    if (dapl < 0)
        return fail(r, "unable to prepare a dataset for reading");

    r->dataset[id] = H5Dopen2(r->file, OUT_FIELDS[id].name, dapl);
    H5Pclose(dapl);

    if (r->dataset[id] < 0)
        return fail_field(r, id, "not present");

    return check_shape(r, id, dims);
}

h5reader *h5reader_open(const char *path)
{
    h5reader *r = calloc(1, sizeof *r);
    if (!r)
        return NULL;

    /* Report failures through h5reader_error rather than HDF5's own stack
     * trace on stderr. */
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++)
        r->dataset[id] = H5I_INVALID_HID;

    r->reads = H5I_INVALID_HID;

    r->file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (r->file < 0) {
        fail(r, "unable to open the file, which may not be HDF5 at all");
        return r;
    }

    if (probe_shape(r) < 0)
        return r;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++)
        if (open_field(r, id) < 0)
            return r;

    r->reads = H5Gopen2(r->file, OUT_READS_GROUP, H5P_DEFAULT);
    if (r->reads < 0)
        fail(r, "the group holding the counts about reads is not present");

    return r;
}

void h5reader_close(h5reader *r)
{
    if (!r)
        return;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++)
        if (r->dataset[id] >= 0)
            H5Dclose(r->dataset[id]);

    if (r->reads >= 0)
        H5Gclose(r->reads);

    if (r->file >= 0)
        H5Fclose(r->file);

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

/* Selects the part of a dataset belonging to one reference: the whole row for
 * a row field, the single value for a scalar one. */
static int select_row(hid_t dataset, int32_t tid, size_t width, int rank,
                      hid_t *filespace, hid_t *memspace)
{
    hsize_t start[OUT_RANK_MAX] = { (hsize_t)tid, 0 };
    hsize_t count[OUT_RANK_MAX] = { 1, (hsize_t)width };
    hsize_t extent              = rank == OUT_RANK_VECTOR ? (hsize_t)width : 1;

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

int h5reader_field(h5reader *r, out_field_id id, int32_t tid, double *values)
{
    size_t width = out_extent(id, r->ref_cap, r->ref_cap);
    hid_t  filespace, memspace;
    herr_t status;

    if (tid < 0 || tid >= r->n_refs)
        return fail(r, "reference index outside the file");

    if (select_row(r->dataset[id], tid, width, out_rank(id),
                   &filespace, &memspace) < 0)
        return fail(r, "unable to select an input row");

    status = H5Dread(r->dataset[id], H5T_NATIVE_DOUBLE, memspace, filespace,
                     H5P_DEFAULT, values);

    H5Sclose(memspace);
    H5Sclose(filespace);

    return status < 0 ? fail(r, "unable to read an input row") : 0;
}

/* ------------------------------------------------------------------------ */
/* Totals                                                                    */
/* ------------------------------------------------------------------------ */

int h5reader_count(h5reader *r, const char *name, size_t *value)
{
    hid_t    dataset = H5Dopen2(r->reads, name, H5P_DEFAULT);
    uint64_t stored  = 0;
    herr_t   status;

    if (dataset < 0) {
        snprintf(r->error, sizeof r->error, "%s/%s: not present",
                 OUT_READS_GROUP, name);
        return -1;
    }

    status = H5Dread(dataset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, &stored);
    H5Dclose(dataset);

    if (status < 0) {
        snprintf(r->error, sizeof r->error, "%s/%s: unable to read it",
                 OUT_READS_GROUP, name);
        return -1;
    }

    *value = (size_t)stored;
    return 0;
}
