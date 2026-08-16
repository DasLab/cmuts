/* h5writer.c -- HDF5 output, one row per reference.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "h5writer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "error.h"
#include "h5layout.h"
#include "version.h"

struct h5writer {
    hid_t   file;
    hid_t   dataset[OUT_N_FIELDS];
    /* A dataspace apiece, kept for the life of the writer. Only the selection
     * differs between one row and the next, and a writer is used from one thread,
     * so the handles are made once and reselected. */
    hid_t   filespace[OUT_N_FIELDS];
    hid_t   memspace;   /* one row of the widest field, selected down to size */
    int32_t n_refs;
    size_t  ref_cap;
    double *padding;    /* NaN, as wide as the longest tail a row can have */
    const bool *wanted; /* the optional fields this run writes; NULL for all of them */
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

/* Returns whether the path holds a file with something in it. An empty file is what
 * mktemp and shell redirection leave behind, and has nothing to lose. */
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
/* What wrote the file                                                       */
/* ------------------------------------------------------------------------ */

/* Writes one string attribute on the root group.
 *
 * A file identifier resolves to the root group, so the attribute is written on the file
 * as a whole. The string is stored at its own length and not padded to a fixed one, so two
 * files written by the same program agree byte for byte. */
static int write_identity(hid_t file, const char *name, const char *value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t type  = H5Tcopy(H5T_C_S1);
    hid_t attr  = H5I_INVALID_HID;
    int   status = -1;

    if (space < 0 || type < 0 || H5Tset_size(type, strlen(value)) < 0) {
        goto done;
    }

    attr = H5Acreate2(file, name, type, space, H5P_DEFAULT, H5P_DEFAULT);

    if (attr >= 0 && H5Awrite(attr, type, value) >= 0) {
        status = 0;
    }

done:
    if (attr >= 0) {
        H5Aclose(attr);
    }
    if (type >= 0) {
        H5Tclose(type);
    }
    if (space >= 0) {
        H5Sclose(space);
    }

    return status;
}


/* Records which program wrote the file, and at which version. output.h names the
 * attributes; the values are what this run holds.
 *
 * Neither value varies between two runs over one input, which is what leaves two such
 * runs identical byte for byte. */
static int stamp_identity(h5writer *w, const char *program)
{
    const char *value[OUT_N_ATTRS] = {
        [OUT_ATTR_PROGRAM] = program,
        [OUT_ATTR_VERSION] = CMUTS_VERSION,
    };

    for (out_attr_id id = 0; id < OUT_N_ATTRS; id++) {
        if (write_identity(w->file, OUT_ATTRIBUTES[id].name, value[id]) < 0) {
            return -1;
        }
    }

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Groups                                                                    */
/* ------------------------------------------------------------------------ */

/* Long enough for any name in OUT_FIELDS. A name outgrowing it fails the create rather
 * than being truncated into a group it does not belong to. */
#define GROUP_PATH_MAX 256

/* Creates the group at this path, or leaves the one already there.
 *
 * Untimed, as the file and the datasets are, or two runs over one input would differ by
 * the moment it was created. The groups are made here for that reason:
 * H5Pset_create_intermediate_group creates them with default properties and so stamps
 * each with a creation time. */
static int create_group(h5writer *w, const char *path)
{
    htri_t found = H5Lexists(w->file, path, H5P_DEFAULT);
    hid_t  gcpl, group;

    if (found < 0) {
        return -1;
    }

    if (found > 0) {
        return 0;
    }

    gcpl = h5layout_untimed_plist(H5P_GROUP_CREATE);
    if (gcpl < 0) {
        return -1;
    }

    group = H5Gcreate2(w->file, path, H5P_DEFAULT, gcpl, H5P_DEFAULT);
    H5Pclose(gcpl);

    if (group < 0) {
        return -1;
    }

    H5Gclose(group);
    return 0;
}

/* Creates the groups a field's name places it in, outermost first, so that each is made
 * inside one already present. */
static int create_field_groups(h5writer *w, const char *name)
{
    char path[GROUP_PATH_MAX];

    for (const char *sep = strchr(name, '/'); sep; sep = strchr(sep + 1, '/')) {
        size_t len = (size_t)(sep - name);

        if (len >= sizeof path) {
            return -1;
        }

        memcpy(path, name, len);
        path[len] = '\0';

        if (create_group(w, path) < 0) {
            return -1;
        }
    }

    return 0;
}

/* Creates every group the field names imply, before any dataset is made in one. The names
 * are where the layout fixes which groups an output holds; nothing else declares them. */
static int create_groups(h5writer *w)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (!out_wanted(id, w->wanted)) {
            continue;
        }

        if (create_field_groups(w, OUT_FIELDS[id].name) < 0) {
            return fail(w, "unable to create a group");
        }
    }

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

    space = rank == 0 ? H5Screate(H5S_SCALAR)
                      : H5Screate_simple(rank, dims, NULL);
    dcpl  = h5layout_creation_plist(id, chunk, rank);
    dapl  = h5layout_access_plist(id, chunk, rank);

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

/* Allocates a writer holding nothing yet, with every handle marked absent.
 *
 * The steps that build the rest may each fail and leave those after them undone, and the
 * writer is closed whatever happened, so it must be safe to close from here onwards: it
 * closes exactly what it opened. Zero, which calloc leaves behind, is a handle HDF5 would
 * accept, hence the marking. */
static h5writer *writer_alloc(int32_t n_refs, size_t ref_cap)
{
    h5writer *w = calloc(1, sizeof *w);

    if (!w) {
        return NULL;
    }

    /* Report failures through h5writer_error, with HDF5's own stack trace on stderr
     * turned off. */
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        w->dataset[id]   = H5I_INVALID_HID;
        w->filespace[id] = H5I_INVALID_HID;
    }

    w->file     = H5I_INVALID_HID;
    w->memspace = H5I_INVALID_HID;

    w->n_refs  = n_refs;
    w->ref_cap = ref_cap;

    return w;
}

/* Fills the row of marks that a reference's unused columns are written from. */
static int build_padding(h5writer *w)
{
    w->padding = calloc(w->ref_cap ? w->ref_cap : 1, sizeof *w->padding);

    if (!w->padding) {
        return fail(w, "out of memory");
    }

    for (size_t i = 0; i < w->ref_cap; i++) {
        w->padding[i] = (double)NAN;
    }

    return 0;
}

/* Prepares the row every write is selected from. */
static int build_memspace(h5writer *w)
{
    w->memspace = h5layout_row_space(w->ref_cap, w->wanted);

    return w->memspace < 0 ? fail(w, "unable to prepare the output") : 0;
}

/* Creates the file, exclusively unless overwrite was requested, so that a file appearing
 * between the check and the create cannot undo the decision. */
static int create_file(h5writer *w, const char *path, bool overwrite)
{
    hid_t fcpl = h5layout_untimed_plist(H5P_FILE_CREATE);

    if (fcpl < 0) {
        return fail(w, "unable to prepare the output file");
    }

    w->file = H5Fcreate(path, overwrite ? H5F_ACC_TRUNC : H5F_ACC_EXCL,
                        fcpl, H5P_DEFAULT);
    H5Pclose(fcpl);

    return w->file < 0 ? fail(w, "unable to create the output file") : 0;
}

/* Creates a dataset per field, and the dataspace each of its rows is selected from. */
static int create_fields(h5writer *w)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (!out_wanted(id, w->wanted)) {
            continue;
        }

        w->dataset[id] = create_field(w, id);
        if (w->dataset[id] < 0) {
            return fail(w, "unable to create a dataset");
        }

        w->filespace[id] = H5Dget_space(w->dataset[id]);
        if (w->filespace[id] < 0) {
            return fail(w, "unable to describe a dataset");
        }
    }

    return 0;
}

h5writer *h5writer_create(const char *path, const char *program, int32_t n_refs,
                          size_t ref_cap, bool overwrite, const bool *wanted)
{
    h5writer *w = writer_alloc(n_refs, ref_cap);

    if (!w) {
        return NULL;
    }

    w->wanted = wanted;

    if (build_padding(w) == 0 &&
        build_memspace(w) == 0 &&
        create_file(w, path, overwrite) == 0 &&
        stamp_identity(w, program) == 0 &&
        create_groups(w) == 0) {
        create_fields(w);
    }

    return w;
}

void h5writer_close(h5writer *w)
{
    if (!w) {
        return;
    }

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (w->filespace[id] >= 0) {
            H5Sclose(w->filespace[id]);
        }
        if (w->dataset[id] >= 0) {
            H5Dclose(w->dataset[id]);
        }
    }

    if (w->memspace >= 0) {
        H5Sclose(w->memspace);
    }

    if (w->file >= 0) {
        H5Fclose(w->file);
    }

    free(w->padding);
    free(w);
}

const char *h5writer_error(const h5writer *w)
{
    return w->error[0] ? w->error : NULL;
}

int h5writer_fail(const h5writer *w, const char *path, char *error, size_t error_len)
{
    const char *why = h5writer_error(w);

    snprintf(error, error_len, "%s: %s", path, why ? why : "unable to write it");
    return -1;
}

/* ------------------------------------------------------------------------ */
/* Rows                                                                      */
/* ------------------------------------------------------------------------ */

static int write_part(h5writer *w, out_field_id id, int32_t tid, size_t from,
                      size_t n, const double *values)
{
    herr_t status;

    if (n == 0) {
        return 0;
    }

    if (h5layout_select_span(w->filespace[id], w->memspace, id, tid, from, n) < 0) {
        return fail(w, "unable to select an output row");
    }

    status = H5Dwrite(w->dataset[id], H5T_NATIVE_DOUBLE, w->memspace,
                      w->filespace[id], H5P_DEFAULT, values);

    return status < 0 ? fail(w, "unable to write an output row") : 0;
}

/* Writes one field's row: the reference's own values, then the mark for the columns
 * past them.
 *
 * The tail is written with the row. The two share a chunk, so marking it now writes to
 * a chunk already open, and returning to it later would mean inflating and deflating
 * that chunk again. A reference as long as
 * the longest has no tail. */
int h5writer_field(h5writer *w, out_field_id id, int32_t tid, size_t len,
                   const double *values)
{
    size_t held  = out_values(id, len, w->ref_cap);
    size_t width = out_values(id, w->ref_cap, w->ref_cap);

    if (tid < 0 || tid >= w->n_refs) {
        return fail(w, "reference index outside the output");
    }

    if (write_part(w, id, tid, 0, held, values) < 0) {
        return -1;
    }

    if (held == width) {
        return 0;
    }

    return write_part(w, id, tid, held, width - held, w->padding);
}

/* Writes the whole of one reference's row, in the type the field is stored as.
 *
 * Every column is given a value, so nothing is left to mark as outside the reference. This
 * is the path for values that were read from a file of the same layout: they are already
 * the stored type, and passing them through the double the accumulator uses would widen
 * and narrow them for nothing. */
int h5writer_row(h5writer *w, out_field_id id, int32_t tid, const void *values)
{
    size_t width = out_values(id, w->ref_cap, w->ref_cap);
    herr_t status;

    if (tid < 0 || tid >= w->n_refs) {
        return fail(w, "reference index outside the output");
    }

    if (h5layout_select_span(w->filespace[id], w->memspace, id, tid, 0, width) < 0) {
        return fail(w, "unable to select an output row");
    }

    status = H5Dwrite(w->dataset[id], h5layout_memory_type(id), w->memspace,
                      w->filespace[id], H5P_DEFAULT, values);

    return status < 0 ? fail(w, "unable to write an output row") : 0;
}

bool h5writer_holds(const h5writer *w, out_field_id id)
{
    return out_wanted(id, w->wanted);
}

int h5writer_block(h5writer *w, out_field_id id, int32_t tid, size_t len,
                   const double *values)
{
    herr_t status;

    if (tid < 0 || tid >= w->n_refs) {
        return fail(w, "reference index outside the output");
    }

    if (h5layout_select_block(w->filespace[id], w->memspace, id, tid, len) < 0) {
        return fail(w, "unable to select an output block");
    }

    status = H5Dwrite(w->dataset[id], H5T_NATIVE_DOUBLE, w->memspace,
                      w->filespace[id], H5P_DEFAULT, values);

    return status < 0 ? fail(w, "unable to write an output block") : 0;
}

/* ------------------------------------------------------------------------ */
/* Totals                                                                    */
/* ------------------------------------------------------------------------ */

/* Writes the whole of a field that belongs to the run and not to a reference.
 *
 * Transferred as an unsigned and not through the double every row passes through, since
 * a run total is counted whole and never accumulated. */
int h5writer_value(h5writer *w, out_field_id id, double value)
{
    float  stored = (float)value;
    herr_t status;

    if (OUT_FIELDS[id].per_ref) {
        return fail(w, "a field with a row per reference holds no single value");
    }

    status = H5Dwrite(w->dataset[id], H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, &stored);

    return status < 0 ? fail(w, "unable to write a value") : 0;
}

int h5writer_total(h5writer *w, out_field_id id, size_t value)
{
    uint64_t stored = value;
    herr_t   status;

    if (OUT_FIELDS[id].per_ref) {
        return fail(w, "a field with a row per reference has no run total");
    }

    status = H5Dwrite(w->dataset[id], H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, &stored);

    return status < 0 ? fail(w, "unable to write a run total") : 0;
}
