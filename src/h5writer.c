/* h5writer.c -- HDF5 output, one row per reference.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "h5writer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>

#include "error.h"
#include "h5layout.h"

struct h5writer {
    hid_t   file;
    hid_t   dataset[OUT_N_FIELDS];
    /* A dataspace apiece, kept for the life of the writer. Only the selection
     * differs between one row and the next, and a writer is used from one thread,
     * so the handles are made once and reselected. */
    hid_t   filespace[OUT_N_FIELDS];
    hid_t   memspace;   /* one row of the widest field, selected down to size */
    hid_t   reads;      /* the group the per-run counts are gathered in */
    int32_t n_refs;
    size_t  ref_cap;
    double *padding;    /* NaN, as wide as the longest tail a row can have */
    char    error[CM_ERROR_MAX];
};

static int fail(h5writer *w, const char *what)
{
    snprintf(w->error, sizeof w->error, "%s", what);
    return -1;
}

/* ------------------------------------------------------------------------ */
/* The path                                                                  */
/* ------------------------------------------------------------------------ */

/* Whether the path holds a file with something in it. An empty file is what mktemp
 * and shell redirection leave behind, and has nothing to lose. */
static bool holds_data(const char *path)
{
    struct stat info;

    return stat(path, &info) == 0 && info.st_size > 0;
}

int h5writer_may_replace(const char *path, bool overwrite, bool *may_replace,
                         char *error, size_t error_len)
{
    if (!overwrite && holds_data(path)) {
        snprintf(error, error_len,
                 "%s already holds data; pass --overwrite to replace it", path);
        return -1;
    }

    /* Where nothing is at the path the create stays exclusive, so a file appearing
     * in between is not quietly replaced. An empty file already there has to be
     * truncated instead. */
    *may_replace = overwrite || access(path, F_OK) == 0;
    return 0;
}

/* ------------------------------------------------------------------------ */
/* Dataset construction                                                      */
/* ------------------------------------------------------------------------ */

static hid_t create_field(h5writer *w, out_field_id id)
{
    hsize_t dims[OUT_RANK_MAX]  = { 0, 0 };
    hsize_t chunk[OUT_RANK_MAX] = { 0, 0 };
    int     rank                = out_rank(id);
    hid_t   space, dcpl, dapl, dataset;

    h5layout_shape(id, w->n_refs, w->ref_cap, dims, chunk);

    space = H5Screate_simple(rank, dims, NULL);
    dcpl  = h5layout_creation_plist(id, chunk, rank);
    dapl  = h5layout_access_plist(chunk, rank);

    if (space < 0 || dcpl < 0 || dapl < 0) {
        dataset = H5I_INVALID_HID;
    } else {
        dataset = H5Dcreate2(w->file, OUT_FIELDS[id].name, h5layout_type(id),
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
    hid_t     fcpl, gcpl;
    h5writer *w = calloc(1, sizeof *w);
    if (!w)
        return NULL;

    /* Report failures through h5writer_error rather than HDF5's own stack trace on
     * stderr. */
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    /* Every handle is marked absent before anything can fail, so that a writer
     * abandoned partway through closes exactly what it opened. Zero, which calloc
     * leaves behind, is a handle HDF5 would accept. */
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        w->dataset[id]   = H5I_INVALID_HID;
        w->filespace[id] = H5I_INVALID_HID;
    }

    w->file     = H5I_INVALID_HID;
    w->reads    = H5I_INVALID_HID;
    w->memspace = H5I_INVALID_HID;

    w->n_refs  = n_refs;
    w->ref_cap = ref_cap;
    w->padding = calloc(ref_cap ? ref_cap : 1, sizeof *w->padding);

    if (!w->padding) {
        fail(w, "out of memory");
        return w;
    }

    for (size_t i = 0; i < ref_cap; i++)
        w->padding[i] = (double)NAN;

    w->memspace = h5layout_row_space(ref_cap);
    if (w->memspace < 0) {
        fail(w, "unable to prepare the output");
        return w;
    }

    fcpl = h5layout_untimed_plist(H5P_FILE_CREATE);
    if (fcpl < 0) {
        fail(w, "unable to prepare the output file");
        return w;
    }

    /* Exclusive unless overwrite was requested, so that a file appearing between
     * the check and the create cannot undo the decision. */
    w->file = H5Fcreate(path, overwrite ? H5F_ACC_TRUNC : H5F_ACC_EXCL,
                        fcpl, H5P_DEFAULT);
    H5Pclose(fcpl);

    if (w->file < 0) {
        fail(w, "unable to create the output file");
        return w;
    }

    /* Created before any dataset inside it, so that a dataset's path names a group
     * already present. Untimed, as the file and the datasets are, or two runs over
     * one input would differ by the moment this was created. */
    gcpl = h5layout_untimed_plist(H5P_GROUP_CREATE);
    if (gcpl < 0) {
        fail(w, "unable to prepare the output");
        return w;
    }

    w->reads = H5Gcreate2(w->file, OUT_READS_GROUP, H5P_DEFAULT, gcpl, H5P_DEFAULT);
    H5Pclose(gcpl);

    if (w->reads < 0) {
        fail(w, "unable to create a group");
        return w;
    }

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        w->dataset[id] = create_field(w, id);
        if (w->dataset[id] < 0) {
            fail(w, "unable to create a dataset");
            return w;
        }

        w->filespace[id] = H5Dget_space(w->dataset[id]);
        if (w->filespace[id] < 0) {
            fail(w, "unable to describe a dataset");
            return w;
        }
    }

    return w;
}

void h5writer_close(h5writer *w)
{
    if (!w)
        return;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (w->filespace[id] >= 0)
            H5Sclose(w->filespace[id]);
        if (w->dataset[id] >= 0)
            H5Dclose(w->dataset[id]);
    }

    if (w->memspace >= 0)
        H5Sclose(w->memspace);

    if (w->reads >= 0)
        H5Gclose(w->reads);

    if (w->file >= 0)
        H5Fclose(w->file);

    free(w->padding);
    free(w);
}

const char *h5writer_error(const h5writer *w)
{
    return w->error[0] ? w->error : NULL;
}

/* ------------------------------------------------------------------------ */
/* Rows                                                                      */
/* ------------------------------------------------------------------------ */

static int write_part(h5writer *w, out_field_id id, int32_t tid, size_t from,
                      size_t n, const double *values)
{
    herr_t status;

    if (n == 0)
        return 0;

    if (h5layout_select_span(w->filespace[id], w->memspace, id, tid, from, n) < 0)
        return fail(w, "unable to select an output row");

    status = H5Dwrite(w->dataset[id], H5T_NATIVE_DOUBLE, w->memspace,
                      w->filespace[id], H5P_DEFAULT, values);

    return status < 0 ? fail(w, "unable to write an output row") : 0;
}

/* Writes one field's row: the reference's own values, then the mark for the columns
 * past them.
 *
 * The tail is written with the row rather than swept up at the end. The two share a
 * chunk, so marking it now writes to a chunk already open where returning to it
 * later would mean inflating and deflating that chunk again. A reference as long as
 * the longest has no tail. */
int h5writer_field(h5writer *w, out_field_id id, int32_t tid, size_t len,
                   const double *values)
{
    size_t extent = out_extent(id, len, w->ref_cap);
    size_t width  = out_extent(id, w->ref_cap, w->ref_cap);

    if (tid < 0 || tid >= w->n_refs)
        return fail(w, "reference index outside the output");

    if (write_part(w, id, tid, 0, extent, values) < 0)
        return -1;

    if (extent == width)
        return 0;

    return write_part(w, id, tid, extent, width - extent, w->padding);
}

/* ------------------------------------------------------------------------ */
/* Totals                                                                    */
/* ------------------------------------------------------------------------ */

/* Writes a run total as a scalar dataset beside the per-reference counts.
 *
 * A dataset rather than an attribute on the group, since h5ls does not show
 * attributes unless asked. Stored as the integer type it is counted in rather than
 * the float its neighbors are narrowed to, being counted whole and never
 * accumulated. */
int h5writer_count(h5writer *w, const char *name, size_t value)
{
    uint64_t stored  = value;
    hid_t    space   = H5Screate(H5S_SCALAR);
    hid_t    dcpl    = h5layout_untimed_plist(H5P_DATASET_CREATE);
    hid_t    dataset = H5I_INVALID_HID;
    herr_t   status  = -1;

    if (space >= 0 && dcpl >= 0)
        dataset = H5Dcreate2(w->reads, name, H5T_STD_U64LE, space,
                             H5P_DEFAULT, dcpl, H5P_DEFAULT);

    if (dataset >= 0)
        status = H5Dwrite(dataset, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, &stored);

    if (dataset >= 0)
        H5Dclose(dataset);
    if (dcpl >= 0)
        H5Pclose(dcpl);
    if (space >= 0)
        H5Sclose(space);

    return status < 0 ? fail(w, "unable to write a run total") : 0;
}
