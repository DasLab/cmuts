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
#include "output.h"

struct h5writer {
    hid_t   file;
    hid_t   dataset[OUT_N_FIELDS];
    hid_t   reads;      /* the group the per-run counts are gathered in */
    int32_t n_refs;
    size_t  ref_cap;
    /* NaN, as wide as the longest tail any row can have, so marking one is a
     * write and not a fill each time. */
    double *padding;
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

/* What is worth refusing is a file with something in it. A path that exists
 * but is empty is what mktemp and shell redirection leave behind, and there is
 * nothing there to lose. */
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

    /* Where nothing is at the path, the create stays exclusive, so a file
     * appearing in between is not quietly replaced. Where an empty one is
     * already there it has to be truncated instead. */
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

    out_shape(id, w->n_refs, w->ref_cap, dims, chunk);

    space = H5Screate_simple(rank, dims, NULL);
    dcpl  = out_layout_plist(id, chunk, rank);
    dapl  = out_access_plist(chunk, rank);

    if (space < 0 || dcpl < 0 || dapl < 0) {
        dataset = H5I_INVALID_HID;
    } else {
        dataset = H5Dcreate2(w->file, OUT_FIELDS[id].name, out_type(id),
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

    /* Report failures through h5writer_error rather than HDF5's own stack
     * trace on stderr. */
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    w->n_refs  = n_refs;
    w->ref_cap = ref_cap;
    w->padding = calloc(ref_cap ? ref_cap : 1, sizeof *w->padding);

    if (!w->padding) {
        fail(w, "out of memory");
        return w;
    }

    for (size_t i = 0; i < ref_cap; i++)
        w->padding[i] = (double)NAN;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++)
        w->dataset[id] = H5I_INVALID_HID;

    w->reads = H5I_INVALID_HID;

    fcpl = out_untimed_plist(H5P_FILE_CREATE);
    if (fcpl < 0) {
        fail(w, "unable to prepare the output file");
        return w;
    }

    /* Exclusive unless overwrite was requested, so that the decision cannot
     * be undone by the file appearing between the check and the create. */
    w->file = H5Fcreate(path, overwrite ? H5F_ACC_TRUNC : H5F_ACC_EXCL,
                        fcpl, H5P_DEFAULT);
    H5Pclose(fcpl);

    if (w->file < 0) {
        fail(w, "unable to create the output file");
        return w;
    }

    /* Named before any dataset inside it, so a path names a group that is
     * already there. Untimed as the file and the datasets are, or two runs over
     * one input would differ by the moment this was made. */
    gcpl = out_untimed_plist(H5P_GROUP_CREATE);
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
    }

    return w;
}

void h5writer_close(h5writer *w)
{
    if (!w)
        return;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++)
        if (w->dataset[id] >= 0)
            H5Dclose(w->dataset[id]);

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

/* Selects the part of a dataset belonging to one reference: the whole row for
 * a scalar field, a span of it for a row field. */
static int select_row(hid_t dataset, int32_t tid, size_t from, size_t width,
                      int rank, hid_t *filespace, hid_t *memspace)
{
    hsize_t start[OUT_RANK_MAX] = { (hsize_t)tid, (hsize_t)from };
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

static int write_part(h5writer *w, out_field_id id, int32_t tid, size_t from,
                      size_t n, const double *values)
{
    hid_t  filespace, memspace;
    herr_t status;

    if (select_row(w->dataset[id], tid, from, n, out_rank(id),
                   &filespace, &memspace) < 0)
        return fail(w, "unable to select an output row");

    status = H5Dwrite(w->dataset[id], H5T_NATIVE_DOUBLE, memspace, filespace,
                      H5P_DEFAULT, values);

    H5Sclose(memspace);
    H5Sclose(filespace);

    return status < 0 ? fail(w, "unable to write an output row") : 0;
}

/* The reference's own values, and then the mark for the columns past them.
 *
 * The tail is written with the row rather than swept up at the end: a row and
 * its tail share a chunk, so marking it now costs a chunk already in hand where
 * returning to it later costs reading, inflating and deflating that chunk
 * again. Where a reference is as long as the longest there is no tail, which on
 * a library of one length is every reference. */
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

/* A run total, as a scalar beside the per-reference counts it belongs with.
 *
 * A dataset and not an attribute on the group: h5ls lists a dataset and passes
 * an attribute over unless asked for one, and a figure nothing shows is a
 * figure nobody finds. It keeps the exact type it is counted in rather than
 * the float its neighbours are narrowed to, being counted whole in the loader
 * and never accumulated. */
int h5writer_count(h5writer *w, const char *name, size_t value)
{
    uint64_t stored  = value;
    hid_t    space   = H5Screate(H5S_SCALAR);
    hid_t    dcpl    = out_untimed_plist(H5P_DATASET_CREATE);
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
