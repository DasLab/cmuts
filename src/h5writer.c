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

#include "chunktally.h"
#include "codec.h"
#include "error.h"
#include "h5layout.h"
#include "version.h"

/* A chunk of one field being gathered. */
typedef struct {
    int64_t        index;   /* -1 while the slot holds nothing */
    unsigned char *values;
} chunk_stage;

struct h5chunk {
    out_field_id    id;
    int64_t         index;
    unsigned char  *values;
    unsigned char  *bytes;      /* where the filtered form goes */
    size_t          raw_bytes;
    size_t          elem;
    size_t          n_bytes;    /* zero until filtered; zero after, the deflate failed */
    struct h5chunk *next;       /* the writer's line of chunks awaiting a taker */
};

/* One field written by whole chunks. Each field has its own chunk shape, since h5layout
 * sizes chunks in bytes, so each keeps its own tally and its own stages. */
typedef struct {
    hsize_t      rows;       /* references per chunk */
    size_t       width;      /* values one row holds */
    size_t       elem;       /* bytes one value occupies */
    size_t       raw_bytes;  /* rows * width * elem */
    chunktally  *tally;      /* set only for a gathered field */
    chunk_stage *stages;
    size_t       n_stages;
} chunkfield;

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
    double *padding;    /* NaN, as wide as the most a row can need */
    const bool *wanted; /* the optional fields this run writes; NULL for all of them */

    chunkfield gathered[OUT_N_FIELDS];
    h5chunk   *ready;        /* finished chunks not yet taken, oldest first */
    h5chunk   *ready_tail;

    char    error[CM_ERROR_MAX];
};

static bool gathers(const h5writer *w, out_field_id id)
{
    return w->gathered[id].tally != NULL;
}

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
/* Chunks in and out                                                         */
/* ------------------------------------------------------------------------ */

static h5chunk *chunk_open(const h5writer *w, out_field_id id, int64_t index,
                           unsigned char *values)
{
    h5chunk *chunk = calloc(1, sizeof *chunk);

    if (!chunk) {
        return NULL;
    }

    chunk->id        = id;
    chunk->index     = index;
    chunk->values    = values;
    chunk->raw_bytes = w->gathered[id].raw_bytes;
    chunk->elem      = w->gathered[id].elem;
    chunk->bytes     = malloc(codec_bound(chunk->raw_bytes));

    if (!chunk->bytes) {
        free(chunk);
        return NULL;
    }

    return chunk;
}

static void chunk_free(h5chunk *chunk)
{
    free(chunk->bytes);
    free(chunk->values);
    free(chunk);
}

static void ready_append(h5writer *w, h5chunk *chunk)
{
    if (w->ready_tail) {
        w->ready_tail->next = chunk;
    } else {
        w->ready = chunk;
    }

    w->ready_tail = chunk;
}

h5chunk *h5writer_take_chunk(h5writer *w)
{
    h5chunk *chunk = w->ready;

    if (!chunk) {
        return NULL;
    }

    w->ready = chunk->next;
    if (!w->ready) {
        w->ready_tail = NULL;
    }

    chunk->next = NULL;
    return chunk;
}

void h5chunk_filter(h5chunk *chunk)
{
    unsigned char *scratch = malloc(codec_scratch(chunk->raw_bytes));

    if (scratch) {
        chunk->n_bytes = codec_encode(chunk->bytes, scratch, chunk->values,
                                      chunk->raw_bytes, chunk->elem);
    }

    free(scratch);
}

int h5writer_write_chunk(h5writer *w, h5chunk *chunk)
{
    hsize_t offset[OUT_RANK_MAX] = { (hsize_t)chunk->index *
                                     w->gathered[chunk->id].rows };
    int     status               = 0;

    if (chunk->n_bytes == 0) {
        status = fail(w, "unable to compress an output chunk");
    } else if (H5Dwrite_chunk(w->dataset[chunk->id], H5P_DEFAULT, 0, offset,
                              chunk->n_bytes, chunk->bytes) < 0) {
        status = fail(w, "unable to write an output chunk");
    }

    chunk_free(chunk);
    return status;
}

/* ------------------------------------------------------------------------ */
/* Gathering rows into chunks                                                */
/* ------------------------------------------------------------------------ */

static chunk_stage *stage_holding(chunkfield *f, int64_t index)
{
    for (size_t i = 0; i < f->n_stages; i++) {
        if (f->stages[i].index == index) {
            return &f->stages[i];
        }
    }

    return NULL;
}

static chunk_stage *stage_spare(chunkfield *f)
{
    for (size_t i = 0; i < f->n_stages; i++) {
        if (f->stages[i].index < 0) {
            return &f->stages[i];
        }
    }

    return NULL;
}

/* How many chunks gather at once is bounded by the pipeline rather than by this file: a
 * chunk stays open only while a reference it holds is still in flight, and the context
 * pool bounds how many those can be. */
static chunk_stage *stage_more(chunkfield *f)
{
    size_t       wanted = f->n_stages ? f->n_stages * 2 : 4;
    chunk_stage *grown  = realloc(f->stages, wanted * sizeof *grown);

    if (!grown) {
        return NULL;
    }

    for (size_t i = f->n_stages; i < wanted; i++) {
        grown[i] = (chunk_stage){ .index = -1 };
    }

    f->stages   = grown;
    f->n_stages = wanted;

    return stage_spare(f);
}

/* Writes the field's fill across a fresh chunk, value by value. */
static void prefill(const chunkfield *f, out_field_id id, unsigned char *values)
{
    out_value fill;

    if (out_fill_value(id, &fill) < 0) {
        return;
    }

    for (size_t i = 0; i < f->raw_bytes; i += f->elem) {
        memcpy(values + i, &fill, f->elem);
    }
}

/* The stage gathering this chunk, begun on its first row. NULL only out of memory. */
static chunk_stage *stage_for(chunkfield *f, out_field_id id, int64_t index)
{
    chunk_stage *s = stage_holding(f, index);

    if (s) {
        return s;
    }

    s = stage_spare(f);
    if (!s) {
        s = stage_more(f);
    }

    if (!s) {
        return NULL;
    }

    s->values = malloc(f->raw_bytes);
    if (!s->values) {
        return NULL;
    }

    prefill(f, id, s->values);
    s->index = index;
    return s;
}

/* The row a reference occupies within its chunk's staged values. */
static unsigned char *stage_row(const chunkfield *f, chunk_stage *s, int32_t tid)
{
    size_t row = (size_t)tid % (size_t)f->rows;

    return s->values + row * f->width * f->elem;
}

/* Narrows accumulated doubles to the type the field is stored as, which the datatype
 * conversion in H5Dwrite would otherwise do on the writing thread. */
static void narrow(void *dst, out_stored stored, const double *src, size_t n)
{
    switch (stored) {
        case OUT_F32: {
            float *to = dst;

            for (size_t i = 0; i < n; i++) {
                to[i] = (float)src[i];
            }
            return;
        }
        case OUT_U64: {
            uint64_t *to = dst;

            for (size_t i = 0; i < n; i++) {
                to[i] = (uint64_t)src[i];
            }
            return;
        }
        case OUT_N_STORED:
            break;
    }
}

/* Pads a staged row, as write_part does with the padding row for a row written directly. */
static void pad_row(unsigned char *row, out_field_id id, size_t from, size_t width)
{
    size_t    elem = out_stored_bytes(id);
    out_value pad;

    if (out_pad_value(id, &pad) < 0) {
        return;
    }

    for (size_t i = from; i < width; i++) {
        memcpy(row + i * elem, &pad, elem);
    }
}

static int64_t gathered_chunk(const chunkfield *f, int32_t tid)
{
    return (int64_t)tid / (int64_t)f->rows;
}

/* Copies one reference's span of a gathered field into its chunk. */
static int gather_span(h5writer *w, out_field_id id, int32_t tid, size_t held,
                       const double *values)
{
    chunkfield  *f = &w->gathered[id];
    chunk_stage *s = stage_for(f, id, gathered_chunk(f, tid));

    if (!s) {
        return fail(w, "out of memory gathering an output chunk");
    }

    narrow(stage_row(f, s, tid), OUT_FIELDS[id].stored, values, held);
    return 0;
}

/* Copies one reference's block of a gathered field into its chunk, a row of the block at
 * a time, at the stride of the full-width row. The columns and rows past the block keep
 * the fill the chunk was prefilled with. */
static int gather_block(h5writer *w, out_field_id id, int32_t tid, size_t len,
                        const double *values)
{
    chunkfield  *f = &w->gathered[id];
    chunk_stage *s = stage_for(f, id, gathered_chunk(f, tid));
    size_t       full[OUT_RANK_MAX];
    size_t       part[OUT_RANK_MAX];
    size_t       stride, rows, cols;
    int          rank;
    unsigned char *row;

    if (!s) {
        return fail(w, "out of memory gathering an output chunk");
    }

    /* A block spans two extents -- there is no other shape wider than one -- so its
     * dataset has the reference dimension and those two, and rank is 3. */
    rank   = out_dims(id, 1, len, part);
    out_dims(id, 1, w->ref_cap, full);
    stride = full[rank - 1];
    cols   = part[rank - 1];
    rows   = part[1];

    row = stage_row(f, s, tid);
    for (size_t i = 0; i < rows; i++) {
        narrow(row + i * stride * f->elem, OUT_FIELDS[id].stored,
               values + i * cols, cols);
    }

    return 0;
}

/* Lines up a finished chunk of one field for h5writer_take_chunk. The buffer travels
 * with the chunk, so the stage is left empty. */
static int send_chunk(h5writer *w, out_field_id id, int64_t index)
{
    chunkfield  *f = &w->gathered[id];
    chunk_stage *s = stage_holding(f, index);
    h5chunk     *chunk;

    if (!s) {
        return 0;
    }

    chunk = chunk_open(w, id, index, s->values);
    if (!chunk) {
        return fail(w, "out of memory sending an output chunk");
    }

    s->values = NULL;
    s->index  = -1;
    ready_append(w, chunk);
    return 0;
}

static int send_settled(h5writer *w)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        int64_t index;

        if (!gathers(w, id)) {
            continue;
        }

        while ((index = chunktally_take_settled(w->gathered[id].tally)) >= 0) {
            if (send_chunk(w, id, index) < 0) {
                return -1;
            }
        }
    }

    return 0;
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

/* Prepares the gathered form of every field whose row holds more than one value. The
 * fields a value wide stay on the row path: their chunks span thousands of references,
 * which would hold nearly the whole output in stages. */
static int build_chunkfields(h5writer *w)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        hsize_t     dims[OUT_RANK_MAX], chunk[OUT_RANK_MAX];
        chunkfield *f     = &w->gathered[id];
        size_t      width = out_values(id, w->ref_cap, w->ref_cap);

        if (!OUT_FIELDS[id].per_ref || !out_wanted(id, w->wanted) || width <= 1) {
            continue;
        }

        h5layout_shape(id, w->n_refs, w->ref_cap, dims, chunk);

        f->rows      = chunk[0];
        f->width     = width;
        f->elem      = out_stored_bytes(id);
        f->raw_bytes = (size_t)f->rows * width * f->elem;
        f->tally     = chunktally_create(w->n_refs, (size_t)f->rows);

        if (!f->tally) {
            return fail(w, "out of memory tracking the output chunks");
        }
    }

    return 0;
}

h5writer *h5writer_create(const char *path, const char *program, int32_t n_refs,
                          size_t ref_cap, bool overwrite, const bool *wanted,
                          bool gather)
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
        create_groups(w) == 0 &&
        create_fields(w) == 0 &&
        gather) {
        build_chunkfields(w);
    }

    return w;
}

/* Finishes every chunk still gathering and writes whatever was never taken, filtering it
 * here, before HDF5 shuts down: a chunk written after the datasets closed would be lost.
 * Chunks a caller took have all come back by now, as the interface requires. */
static void finish_chunks(h5writer *w)
{
    h5chunk *chunk;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (gathers(w, id)) {
            chunktally_no_more(w->gathered[id].tally);
        }
    }

    send_settled(w);

    while ((chunk = h5writer_take_chunk(w))) {
        h5chunk_filter(chunk);
        h5writer_write_chunk(w, chunk);
    }

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        chunkfield *f = &w->gathered[id];

        for (size_t i = 0; i < f->n_stages; i++) {
            free(f->stages[i].values);
        }

        free(f->stages);
        chunktally_destroy(f->tally);
    }
}

void h5writer_close(h5writer *w)
{
    if (!w) {
        return;
    }

    finish_chunks(w);

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

int h5writer_field(h5writer *w, out_field_id id, int32_t tid, size_t len,
                   const double *values)
{
    size_t held = out_values(id, len, w->ref_cap);

    if (tid < 0 || tid >= w->n_refs) {
        return fail(w, "reference index outside the output");
    }

    if (gathers(w, id)) {
        return gather_span(w, id, tid, held, values);
    }

    return write_part(w, id, tid, 0, held, values);
}

int h5writer_pad(h5writer *w, out_field_id id, int32_t tid, size_t len)
{
    size_t held  = out_values(id, len, w->ref_cap);
    size_t width = out_values(id, w->ref_cap, w->ref_cap);

    if (tid < 0 || tid >= w->n_refs) {
        return fail(w, "reference index outside the output");
    }

    if (held == width) {
        return 0;
    }

    if (gathers(w, id)) {
        chunkfield  *f = &w->gathered[id];
        chunk_stage *s = stage_for(f, id, gathered_chunk(f, tid));

        if (!s) {
            return fail(w, "out of memory gathering an output chunk");
        }

        pad_row(stage_row(f, s, tid), id, held, f->width);
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

void h5writer_expect(h5writer *w, int32_t tid)
{
    if (tid < 0 || tid >= w->n_refs) {
        return;
    }

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (gathers(w, id)) {
            chunktally_expect(w->gathered[id].tally, tid);
        }
    }
}

int h5writer_wrote(h5writer *w, int32_t tid)
{
    if (tid < 0 || tid >= w->n_refs) {
        return 0;
    }

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (gathers(w, id)) {
            chunktally_wrote(w->gathered[id].tally, tid);
        }
    }

    return send_settled(w);
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

    if (gathers(w, id)) {
        return gather_block(w, id, tid, len, values);
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

/* Writes the whole of a field that belongs to the run and not to a reference. */
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

/* Transferred as an unsigned and not through the double every row passes through, since
 * a run total is counted whole and never accumulated. */
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
