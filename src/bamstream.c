/* bamstream.c -- a merge across the files, keyed on the reference.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "bamstream.h"

#include <stdio.h>
#include <stdlib.h>

#include <htslib/thread_pool.h>

#include "error.h"

/* Where a reference stands in the order the files are drained. Unmapped reads sit at
 * the end of a coordinate-sorted file, so they sort after every reference rather
 * than before them all, as their tid of -1 would otherwise put them. */
#define ORDER_UNMAPPED INT64_MAX
#define ORDER_NONE     (-1)

/* One file, and the record read ahead from it. The record itself stays in the reader
 * that produced it, so the merge needs no storage of its own. */
typedef struct {
    cm_bam_reader *reader;
    const char    *path;     /* borrowed from the command line */
    cm_bam_record  pending;
    bool           spent;
} source;

struct cm_bam_stream {
    source       *sources;
    size_t        n;
    size_t        at;        /* the source now supplying records */
    int64_t       current;   /* the reference being drained */
    bool          supplied;  /* whether that source's record has been handed out */
    htsThreadPool pool;
    char          error[CM_ERROR_MAX];
};

/* A failure reports the file at fault, which is not otherwise clear. */
static int fail(cm_bam_stream *stream, const char *path, const char *what)
{
    snprintf(stream->error, sizeof stream->error, "%s: %s", path, what);
    return -1;
}

/* ------------------------------------------------------------------------ */
/* Lifetime                                                                  */
/* ------------------------------------------------------------------------ */

static int open_source(cm_bam_stream *stream, source *src, const char *path,
                       const char *fasta_path)
{
    const char *why;

    src->path   = path;
    src->reader = cm_bam_open(path, &why);

    if (!src->reader) {
        return fail(stream, path, why);
    }

    if (!cm_bam_is_coordinate_sorted(src->reader)) {
        return fail(stream, path, "not coordinate sorted; references would stay "
                                  "live to the end of the file");
    }

    if (cm_bam_nref(src->reader) < 1) {
        return fail(stream, path, "declares no references");
    }

    /* Set before any record is read, so that a CRAM decodes against the same
     * reference its reads are compared to. */
    if (cm_bam_set_reference(src->reader, fasta_path) < 0) {
        return fail(stream, path, cm_bam_error(src->reader));
    }

    return 0;
}

/* Points every reader at one shared thread pool, so that the threads follow whichever
 * file is being drained rather than each file holding a set that idles while the
 * others are read. */
static int share_threads(cm_bam_stream *stream, int threads)
{
    if (threads < 1) {
        return 0;
    }

    stream->pool.pool = hts_tpool_init(threads);
    if (!stream->pool.pool) {
        snprintf(stream->error, sizeof stream->error,
                 "unable to start the decompression threads");
        return -1;
    }

    for (size_t i = 0; i < stream->n; i++) {
        if (cm_bam_use_pool(stream->sources[i].reader, &stream->pool) < 0) {
            return fail(stream, stream->sources[i].path,
                        cm_bam_error(stream->sources[i].reader));
        }
    }

    return 0;
}

/* Reads the record a source will offer next, marking it spent at end of file. A spent
 * source takes no further part in the merge. */
static int advance(cm_bam_stream *stream, source *src)
{
    int status = cm_bam_next(src->reader, &src->pending);

    if (status == CM_ITER_ERROR) {
        return fail(stream, src->path, cm_bam_error(src->reader));
    }

    if (status == CM_ITER_EOF) {
        src->spent = true;
    }

    return 0;
}

cm_bam_stream *cm_bam_stream_open(const char *const *paths, size_t n_paths,
                                  const char *fasta_path, int threads)
{
    cm_bam_stream *stream = calloc(1, sizeof *stream);
    if (!stream) {
        return NULL;
    }

    stream->sources = calloc(n_paths, sizeof *stream->sources);
    if (!stream->sources) {
        free(stream);
        return NULL;
    }

    stream->n       = n_paths;
    stream->current = ORDER_NONE;

    for (size_t i = 0; i < n_paths; i++) {
        if (open_source(stream, &stream->sources[i], paths[i], fasta_path) < 0) {
            return stream;
        }
    }

    if (share_threads(stream, threads) < 0) {
        return stream;
    }

    for (size_t i = 0; i < n_paths; i++) {
        if (advance(stream, &stream->sources[i]) < 0) {
            return stream;
        }
    }

    return stream;
}

void cm_bam_stream_close(cm_bam_stream *stream)
{
    if (!stream) {
        return;
    }

    for (size_t i = 0; i < stream->n; i++) {
        cm_bam_close(stream->sources[i].reader);
    }

    /* Destroyed after the readers, which draw on it until they are closed. */
    if (stream->pool.pool) {
        hts_tpool_destroy(stream->pool.pool);
    }

    free(stream->sources);
    free(stream);
}

const char *cm_bam_stream_error(const cm_bam_stream *stream)
{
    return stream->error[0] ? stream->error : NULL;
}

/* ------------------------------------------------------------------------ */
/* The merge                                                                 */
/* ------------------------------------------------------------------------ */

static int64_t order_of(const source *src)
{
    return src->pending.tid == CM_BAM_NO_REF ? ORDER_UNMAPPED : src->pending.tid;
}

/* The first source at or after `from` still holding the given reference. */
static bool holder_of(const cm_bam_stream *stream, int64_t reference, size_t from,
                      size_t *out)
{
    for (size_t i = from; i < stream->n; i++) {
        if (!stream->sources[i].spent && order_of(&stream->sources[i]) == reference) {
            *out = i;
            return true;
        }
    }

    return false;
}

/* The source holding the lowest reference any of them has left. Ties go to the first,
 * so the files are drained in the order given on the command line. */
static bool lowest_holder(const cm_bam_stream *stream, size_t *out)
{
    size_t lowest = 0;
    bool   found  = false;

    for (size_t i = 0; i < stream->n; i++) {
        if (stream->sources[i].spent) {
            continue;
        }

        if (!found || order_of(&stream->sources[i]) < order_of(&stream->sources[lowest])) {
            lowest = i;
            found  = true;
        }
    }

    if (found) {
        *out = lowest;
    }

    return found;
}

/* Chooses the source the next record comes from: the current one while it still holds
 * the reference being drained, then the next source holding it, and only once none
 * does the lowest reference left anywhere. */
static bool select_source(cm_bam_stream *stream)
{
    size_t at;

    if (holder_of(stream, stream->current, stream->at, &at)) {
        stream->at = at;
        return true;
    }

    if (!lowest_holder(stream, &at)) {
        return false;
    }

    stream->at      = at;
    stream->current = order_of(&stream->sources[at]);

    return true;
}

int cm_bam_stream_next(cm_bam_stream *stream, cm_bam_record *out)
{
    /* The record handed out last time belongs to the source, so it is overwritten only
     * now that the caller is done with it. */
    if (stream->supplied && advance(stream, &stream->sources[stream->at]) < 0) {
        return CM_ITER_ERROR;
    }

    stream->supplied = false;

    if (!select_source(stream)) {
        return CM_ITER_EOF;
    }

    stream->supplied = true;
    *out             = stream->sources[stream->at].pending;

    return CM_ITER_OK;
}

const bam1_t *cm_bam_stream_raw(const cm_bam_stream *stream)
{
    return cm_bam_raw(stream->sources[stream->at].reader);
}

/* ------------------------------------------------------------------------ */
/* Position                                                                  */
/* ------------------------------------------------------------------------ */

uint64_t cm_bam_stream_position(const cm_bam_stream *stream)
{
    uint64_t total = 0;

    for (size_t i = 0; i < stream->n; i++) {
        total += cm_bam_position(stream->sources[i].reader);
    }

    return total;
}

/* The summed span, or 0 where the size of any one file is unknown. */
uint64_t cm_bam_stream_span(const cm_bam_stream *stream)
{
    uint64_t total = 0;

    for (size_t i = 0; i < stream->n; i++) {
        uint64_t span = cm_bam_span(stream->sources[i].reader);

        if (span == 0) {
            return 0;
        }

        total += span;
    }

    return total;
}

/* ------------------------------------------------------------------------ */
/* Header queries                                                            */
/* ------------------------------------------------------------------------ */

static const cm_bam_reader *first_header(const cm_bam_stream *stream)
{
    return stream->sources[0].reader;
}

const char *cm_bam_stream_refname(const cm_bam_stream *stream, int32_t tid)
{
    return cm_bam_refname(first_header(stream), tid);
}

int32_t cm_bam_stream_nref(const cm_bam_stream *stream)
{
    return cm_bam_nref(first_header(stream));
}

hts_pos_t cm_bam_stream_reflen(const cm_bam_stream *stream, int32_t tid)
{
    return cm_bam_reflen(first_header(stream), tid);
}

hts_pos_t cm_bam_stream_max_reflen(const cm_bam_stream *stream)
{
    return cm_bam_max_reflen(first_header(stream));
}

size_t cm_bam_stream_count(const cm_bam_stream *stream)
{
    return stream->n;
}

const cm_bam_reader *cm_bam_stream_reader(const cm_bam_stream *stream, size_t file)
{
    return stream->sources[file].reader;
}

const char *cm_bam_stream_path(const cm_bam_stream *stream, size_t file)
{
    return stream->sources[file].path;
}
