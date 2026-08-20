/* h5reader.c -- HDF5 input, one row per reference.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "h5reader.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "h5layout.h"

/* The field whose shape gives the rest. Reactivity is one value per base and required by
 * every program that reads an output, so the width of its row is the capacity every
 * other field's width derives from. */
#define SHAPE_FIELD FMT_REACTIVITY

struct h5reader {
    hid_t   file;
    hid_t   dataset[FMT_N_FIELDS];
    /* A dataspace apiece, kept for the life of the reader as the writer keeps its
     * own: only the selection differs between one row and the next. */
    hid_t   filespace[FMT_N_FIELDS];
    bool    taken[FMT_N_FIELDS];     /* the fields it opened */
    bool    required[FMT_N_FIELDS];  /* the fields the caller's requests insist on */
    hid_t   memspace;   /* one row of the widest field, selected down to size */
    int32_t n_refs;
    size_t  ref_cap;
    char    ignored[CM_ERROR_MAX];  /* datasets the file holds and open_fields did not open */
    char    error[CM_ERROR_MAX];
};

static int fail(h5reader *r, const char *what)
{
    snprintf(r->error, sizeof r->error, "%s", what);
    return -1;
}

static int fail_field(h5reader *r, fmt_field_id id, const char *what)
{
    snprintf(r->error, sizeof r->error, "%s: %s", FMT_FIELDS[id].name, what);
    return -1;
}

/* ------------------------------------------------------------------------ */
/* Shape                                                                     */
/* ------------------------------------------------------------------------ */

/* Returns a dataset's dimensions, failing where its rank is not the one expected. */
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
    hid_t   dataset = H5Dopen2(r->file, FMT_FIELDS[SHAPE_FIELD].name, H5P_DEFAULT);
    hsize_t dims[FMT_RANK_MAX];
    int     status;

    if (dataset < 0) {
        return fail_field(r, SHAPE_FIELD, "not present; this is not a cmuts output");
    }

    status = dataset_dims(dataset, fmt_rank(SHAPE_FIELD), dims);
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

static int check_shape(h5reader *r, fmt_field_id id, const hsize_t *expected)
{
    int     rank = fmt_rank(id);
    hsize_t dims[FMT_RANK_MAX];

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

static int open_field(h5reader *r, fmt_field_id id)
{
    hsize_t dims[FMT_RANK_MAX]  = { 0, 0 };
    hsize_t chunk[FMT_RANK_MAX] = { 0, 0 };
    hid_t   dapl;

    h5layout_shape(id, r->n_refs, r->ref_cap, dims, chunk);

    dapl = h5layout_access_plist(id, chunk, fmt_rank(id));
    if (dapl < 0) {
        return fail(r, "unable to prepare a dataset for reading");
    }

    r->dataset[id] = H5Dopen2(r->file, FMT_FIELDS[id].name, dapl);
    H5Pclose(dapl);

    /* A field no request requires is skipped where the file lacks it. */
    if (r->dataset[id] < 0) {
        if (r->required[id]) {
            return fail_field(r, id, "not present");
        }

        r->taken[id] = false;
        return 0;
    }

    if (check_shape(r, id, dims) < 0) {
        return -1;
    }

    r->filespace[id] = H5Dget_space(r->dataset[id]);

    return r->filespace[id] < 0 ? fail_field(r, id, "cannot be described") : 0;
}

/* Opens the requested fields. */
static int open_fields(h5reader *r)
{
    for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
        if (!fmt_wanted(id, r->taken)) {
            continue;
        }

        if (open_field(r, id) < 0) {
            return -1;
        }
    }

    return 0;
}

/* Appends one dataset path to the list of ignored datasets, comma separated. */
static herr_t note_ignored(hid_t obj, const char *name, const H5O_info2_t *info,
                           void *op_data)
{
    h5reader *r    = op_data;
    size_t    used = strlen(r->ignored);

    (void)obj;

    if (info->type != H5O_TYPE_DATASET) {
        return 0;
    }

    for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
        if (r->taken[id] && strcmp(name, FMT_FIELDS[id].name) == 0) {
            return 0;
        }
    }

    snprintf(r->ignored + used, sizeof r->ignored - used, "%s%s",
             used ? ", " : "", name);
    return 0;
}

/* Collects the datasets the reader did not open. Best effort: a walk that fails leaves
 * the list short, and the fields themselves have been checked already. */
static void find_ignored(h5reader *r)
{
    H5Ovisit3(r->file, H5_INDEX_NAME, H5_ITER_NATIVE, note_ignored, r,
              H5O_INFO_BASIC);
}

/* Allocates a reader holding no handles yet, every one marked absent.
 *
 * The steps that build the rest may each fail and leave those after them undone, and the
 * reader is always closed, so it must be safe to close from here onwards: it
 * closes exactly what it opened. Zero, which calloc leaves behind, is a handle HDF5 would
 * accept, hence the marking. */
static h5reader *reader_alloc(const fmt_request *requests, size_t n)
{
    h5reader *r = calloc(1, sizeof *r);

    if (!r) {
        return NULL;
    }

    for (size_t i = 0; i < n; i++) {
        r->taken[requests[i].id]    = true;
        r->required[requests[i].id] = requests[i].required;
    }

    /* Report failures through h5reader_error, with HDF5's own stack trace on stderr
     * turned off. */
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
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
    r->memspace = h5layout_row_space(r->ref_cap, r->taken);

    return r->memspace < 0 ? fail(r, "unable to prepare the file for reading") : 0;
}

h5reader *h5reader_open(const char *path, const fmt_request *requests, size_t n)
{
    h5reader *r = reader_alloc(requests, n);

    if (!r) {
        return NULL;
    }

    if (open_file(r, path) == 0 &&
        probe_shape(r) == 0 &&
        build_memspace(r) == 0 &&
        open_fields(r) == 0) {
        find_ignored(r);
    }

    return r;
}

void h5reader_close(h5reader *r)
{
    if (!r) {
        return;
    }

    for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
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

const char *h5reader_ignored(const h5reader *r)
{
    return r->ignored[0] ? r->ignored : NULL;
}

int h5reader_fail(const h5reader *r, const char *path, char *error, size_t error_len)
{
    const char *why = h5reader_error(r);

    snprintf(error, error_len, "%s: %s", path, why ? why : "unable to read it");
    return -1;
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

bool h5reader_holds(const h5reader *r, fmt_field_id id)
{
    (void)r;

    return fmt_wanted(id, r->taken);
}

int h5reader_field(h5reader *r, fmt_field_id id, int32_t tid, void *values)
{
    size_t width = fmt_values(id, r->ref_cap, r->ref_cap);
    herr_t status;

    if (tid < 0 || tid >= r->n_refs) {
        return fail(r, "reference index outside the file");
    }

    if (h5layout_select_span(r->filespace[id], r->memspace, id, tid, width) < 0) {
        return fail(r, "unable to select an input row");
    }

    status = H5Dread(r->dataset[id], h5layout_memory_type(id), r->memspace,
                     r->filespace[id], H5P_DEFAULT, values);

    return status < 0 ? fail(r, "unable to read an input row") : 0;
}

/* ------------------------------------------------------------------------ */
/* Totals                                                                    */
/* ------------------------------------------------------------------------ */

int h5reader_total(h5reader *r, fmt_field_id id, size_t *value)
{
    uint64_t stored = 0;
    herr_t   status;

    if (FMT_FIELDS[id].per_ref) {
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
