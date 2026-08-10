/* h5writer.c -- HDF5 output, one row per reference.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "h5writer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <hdf5.h>

#include "error.h"
#include "rates.h"

#define DATASET_REFERENCE "reference"

/* Chunks are sized by bytes rather than rows so that a file of few long
 * references and one of many short references both land near this figure. */
#define TARGET_CHUNK_BYTES (1u << 20)

/* Rows are completed in roughly ascending order, so holding a handful of
 * chunks in the cache keeps a single-row write from forcing a read, modify and
 * rewrite of the chunk around it. The slots are a hash table over chunk
 * indices, which is why their number is prime. */
#define CACHED_CHUNKS 4
#define CACHE_SLOTS   521

/* HDF5's own default, restated because naming either of the other two means
 * passing all three. */
#define CACHE_PREEMPTION 0.75

#define DEFLATE_LEVEL 3

/* A scalar field is one value per reference; every other kind is a reference by
 * something, and so a row. Both shapes are described through the same arrays,
 * which are sized to the larger. */
/* What the output holds, which is no longer what the accumulator holds. The
 * mutations and the span are the evidence gathered; what a caller wants is the
 * rate they come to and how far it is to be believed, so those are what is
 * written and the two they are made of stay inside. */
typedef enum {
    OUT_COVERAGE,
    OUT_REACTIVITY,
    OUT_ERROR,
    OUT_LENGTHS,
    OUT_READS,
    OUT_REJECTED,
    OUT_N_FIELDS,
} out_field_id;

/* A count is whole however it was accumulated: the arena is doubles throughout,
 * so that any two accumulators merge, and what is only ever a tally of ones
 * comes back out of it as one. Written as a float it would round above two to
 * the twenty-fourth, which a deeply read reference passes. What is fractional
 * in its own right -- the coverage, the rate, its error -- stays float. */
typedef struct {
    const char    *name;
    accum_field_id shape;   /* the accumulated field whose extent it takes */
    bool           counted; /* whole, and written as an unsigned */
} out_field;

static const out_field OUT_FIELDS[OUT_N_FIELDS] = {
    [OUT_COVERAGE]   = { "coverage",       ACCUM_COVERAGE, false },
    [OUT_REACTIVITY] = { "reactivity",     ACCUM_MUTATIONS, false },
    [OUT_ERROR]      = { "error",          ACCUM_MUTATIONS, false },
    [OUT_LENGTHS]    = { "reads/lengths",  ACCUM_LENGTHS, true },
    [OUT_READS]      = { "reads/counted",  ACCUM_READS, true },
    [OUT_REJECTED]   = { "reads/rejected", ACCUM_FILTERED, true },
};

#define RANK_SCALAR 1
#define RANK_VECTOR 2
#define RANK_MAX    RANK_VECTOR

struct h5writer {
    hid_t   file;
    hid_t   dataset[OUT_N_FIELDS];
    int32_t n_refs;
    size_t  ref_cap;
    /* NaN, as wide as the longest tail any row can have, so marking one is a
     * write and not a fill each time. */
    hid_t   reads;      /* the group the per-run counts are gathered in */
    double *padding;
    double     *row;    /* what a derived field is worked out into */
    rate_config rates;
    char    error[CM_ERROR_MAX];
};

static int fail(h5writer *w, const char *what)
{
    snprintf(w->error, sizeof w->error, "%s", what);
    return -1;
}

/* A creation property list with object timestamping turned off.
 *
 * HDF5 stamps every object header with the time it was written, so two runs
 * over the same input would produce files differing in bytes that say nothing
 * about the result. Applied to the file as well as the datasets: an fcpl also
 * carries the root group's creation properties, and that group is stamped like
 * any other. */
static hid_t untimed_plist(hid_t class_id)
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

/* ------------------------------------------------------------------------ */
/* Dataset construction                                                      */
/* ------------------------------------------------------------------------ */

static hsize_t rows_per_chunk(const h5writer *w, size_t row_values)
{
    size_t  row_bytes = row_values * sizeof(float);
    hsize_t rows      = row_bytes ? TARGET_CHUNK_BYTES / row_bytes : (hsize_t)w->n_refs;

    if (rows < 1)
        rows = 1;
    if (rows > (hsize_t)w->n_refs)
        rows = (hsize_t)w->n_refs;

    return rows;
}

static size_t out_extent(out_field_id id, size_t len, size_t cap)
{
    return accum_extent(OUT_FIELDS[id].shape, len, cap);
}

static int field_rank(out_field_id id)
{
    return ACCUM_FIELDS[OUT_FIELDS[id].shape].kind == ACCUM_SCALAR
         ? RANK_SCALAR : RANK_VECTOR;
}

/* The width of a row is the field's own extent at the longest reference, so a
 * field wider than one value per base is sized by the same rule as the rest. */
static void field_shape(const h5writer *w, out_field_id id, hsize_t *dims, hsize_t *chunk)
{
    size_t width = out_extent(id, w->ref_cap, w->ref_cap);

    dims[0]  = (hsize_t)w->n_refs;
    chunk[0] = rows_per_chunk(w, width);

    if (field_rank(id) == RANK_VECTOR) {
        dims[1]  = width;
        chunk[1] = width;
    }
}

/* Zero for what is counted, and NaN for what is derived from it.
 *
 * A position no read reached was reached by no read, which is what a count of
 * zero says. A rate is not a count: a reference no read named has no rate, and
 * filling one with zero would say its every position was measured and found
 * unmodified, which is the most confident thing the output can say and it would
 * be saying it about nothing at all. */
static hid_t field_type(out_field_id id)
{
    return OUT_FIELDS[id].counted ? H5T_STD_U64LE : H5T_IEEE_F32LE;
}

/* Zero for a count and for what is measured, NaN for a rate that was not. A
 * count has no NaN to be had, being an unsigned, and needs none: nothing it is
 * written for has padding. */
static const void *field_fill(out_field_id id)
{
    static const uint64_t none = 0;
    static const float    zero = 0.0f;
    static const float    nan  = (float)NAN;

    if (OUT_FIELDS[id].counted)
        return &none;

    return id == OUT_REACTIVITY || id == OUT_ERROR ? (const void *)&nan
                                                   : (const void *)&zero;
}

/* Zero, for every counted field. A position no read reached was reached by no read,
 * which is what a count of zero says, and a reference no read named is only
 * that case for all of its bases at once. What zero must not be taken for is a
 * position outside the reference altogether; those are marked NaN as the row is
 * written, there being no second fill value to say it. */
#define FILL_VALUE 0.0f

static hid_t make_layout(const hsize_t *chunk, int rank, hid_t type,
                         const void *fill)
{
    hid_t dcpl = untimed_plist(H5P_DATASET_CREATE);

    if (dcpl < 0)
        return H5I_INVALID_HID;

    if (H5Pset_chunk(dcpl, rank, chunk) < 0 ||
        H5Pset_fill_value(dcpl, type, fill) < 0 ||
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

    if (H5Pset_chunk_cache(dapl, CACHE_SLOTS, CACHED_CHUNKS * bytes,
                           CACHE_PREEMPTION) < 0) {
        H5Pclose(dapl);
        return H5I_INVALID_HID;
    }

    return dapl;
}

static hid_t create_field(h5writer *w, out_field_id id)
{
    hsize_t dims[RANK_MAX]  = { 0, 0 };
    hsize_t chunk[RANK_MAX] = { 0, 0 };
    int     rank            = field_rank(id);
    hid_t   space, dcpl, dapl, dataset;

    field_shape(w, id, dims, chunk);

    space = H5Screate_simple(rank, dims, NULL);
    dcpl  = make_layout(chunk, rank, field_type(id), field_fill(id));
    dapl  = make_access(chunk, rank);

    if (space < 0 || dcpl < 0 || dapl < 0) {
        dataset = H5I_INVALID_HID;
    } else {
        dataset = H5Dcreate2(w->file, OUT_FIELDS[id].name, field_type(id),
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
                          rate_config rates, bool overwrite)
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
    w->rates     = rates;
    w->padding   = calloc(ref_cap ? ref_cap : 1, sizeof *w->padding);
    w->row       = calloc(ref_cap ? ref_cap : 1, sizeof *w->row);

    if (!w->padding || !w->row) {
        fail(w, "out of memory");
        return w;
    }

    for (size_t i = 0; i < ref_cap; i++)
        w->padding[i] = (double)NAN;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++)
        w->dataset[id] = H5I_INVALID_HID;

    w->reads = H5I_INVALID_HID;

    fcpl = untimed_plist(H5P_FILE_CREATE);
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
    gcpl = untimed_plist(H5P_GROUP_CREATE);
    if (gcpl < 0) {
        fail(w, "unable to prepare the output");
        return w;
    }

    w->reads = H5Gcreate2(w->file, "reads", H5P_DEFAULT, gcpl, H5P_DEFAULT);
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

    free(w->row);
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
 * a scalar field, its first len values for a per-base field. */
static int select_row(hid_t dataset, int32_t tid, size_t from, size_t width,
                      int rank, hid_t *filespace, hid_t *memspace)
{
    hsize_t start[RANK_MAX] = { (hsize_t)tid, (hsize_t)from };
    hsize_t count[RANK_MAX] = { 1, (hsize_t)width };
    hsize_t extent          = rank == RANK_VECTOR ? (hsize_t)width : 1;

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

    if (select_row(w->dataset[id], tid, from, n, field_rank(id),
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
/* What a field is written from. Two of them are not accumulated at all: the
 * rate and its error are what the mutations and the span come to, which rates
 * works out into a row of its own. */
static const double *values(h5writer *w, out_field_id id, const accum *acc,
                            size_t len)
{
    if (id == OUT_REACTIVITY)
        rate_reactivity(&w->rates, acc, len, w->row);
    else if (id == OUT_ERROR)
        rate_error(&w->rates, acc, len, w->row);
    else
        return accum_const_data(acc, OUT_FIELDS[id].shape);

    return w->row;
}

/* The reference's own values, and then the mark for the columns past them.
 *
 * The tail is written with the row rather than swept up at the end: a row and
 * its tail share a chunk, so marking it now costs a chunk already in hand where
 * returning to it later costs reading, inflating and deflating that chunk
 * again. Where a reference is as long as the longest there is no tail, which on
 * a library of one length is every reference. */
static int write_field(h5writer *w, out_field_id id, int32_t tid, size_t len,
                       const accum *acc)
{
    size_t extent = out_extent(id, len, w->ref_cap);
    size_t width  = out_extent(id, w->ref_cap, w->ref_cap);

    if (write_part(w, id, tid, 0, extent, values(w, id, acc, extent)) < 0)
        return -1;

    if (extent == width)
        return 0;

    return write_part(w, id, tid, extent, width - extent, w->padding);
}

int h5writer_row(h5writer *w, int32_t tid, size_t len, const accum *acc)
{
    if (tid < 0 || tid >= w->n_refs)
        return fail(w, "reference index outside the output");

    if (len == 0)
        return 0;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++)
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
    hid_t   dcpl    = untimed_plist(H5P_DATASET_CREATE);
    hid_t   dataset = H5I_INVALID_HID;
    herr_t  status  = -1;

    if (type >= 0 && space >= 0 && dcpl >= 0)
        dataset = H5Dcreate2(w->file, DATASET_REFERENCE, type, space,
                             H5P_DEFAULT, dcpl, H5P_DEFAULT);

    if (dataset >= 0)
        status = H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, names);

    if (dataset >= 0)
        H5Dclose(dataset);
    if (dcpl >= 0)
        H5Pclose(dcpl);
    if (space >= 0)
        H5Sclose(space);
    if (type >= 0)
        H5Tclose(type);

    return status < 0 ? fail(w, "unable to write the reference names") : 0;
}

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
    hid_t    dcpl    = untimed_plist(H5P_DATASET_CREATE);
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
