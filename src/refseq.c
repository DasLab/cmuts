/* refseq.c -- lockstep FASTA streaming, validated against the BAM header.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "refseq.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "checksum.h"
#include "error.h"

struct refseq_source {
    cm_fasta_reader     *fasta;
    const cm_bam_stream *bams;
    cm_fasta_record      record;        /* the record for ordinal loaded - 1 */
    int32_t              loaded;        /* number of records consumed so far */
    cm_bam_sq_cursor    *declarations;  /* each file's @SQ lines, in step */
    char                 error[CM_ERROR_MAX];
};

static const cm_bam_reader *header_of(const refseq_source *src, size_t file)
{
    return cm_bam_stream_reader(src->bams, file);
}

static const char *path_of(const refseq_source *src, size_t file)
{
    return cm_bam_stream_path(src->bams, file);
}

refseq_source *refseq_open(const char *fasta_path, const cm_bam_stream *bams,
                           const char **why)
{
    size_t         files = cm_bam_stream_count(bams);
    refseq_source *src   = calloc(1, sizeof *src);

    if (!src) {
        *why = "out of memory";
        return NULL;
    }

    src->bams  = bams;
    src->fasta = cm_fasta_open(fasta_path, why);

    if (!src->fasta) {
        refseq_close(src);
        return NULL;
    }

    src->declarations = calloc(files, sizeof *src->declarations);

    if (!src->declarations) {
        *why = "out of memory";
        refseq_close(src);
        return NULL;
    }

    for (size_t file = 0; file < files; file++)
        cm_bam_sq_open(&src->declarations[file], header_of(src, file));

    return src;
}

void refseq_close(refseq_source *src)
{
    if (!src)
        return;

    cm_fasta_close(src->fasta);
    free(src->declarations);
    free(src);
}

/* ------------------------------------------------------------------------ */
/* Validation                                                                */
/* ------------------------------------------------------------------------ */

static bool name_matches(refseq_source *src, size_t file, int32_t tid)
{
    const char *expected = cm_bam_refname(header_of(src, file), tid);

    if (expected && strcmp(src->record.name, expected) == 0)
        return true;

    snprintf(src->error, sizeof src->error,
             "reference %d is \"%s\" in the header of %s but \"%s\" in the FASTA; "
             "the two must be in the same order",
             tid, expected ? expected : "(absent)", path_of(src, file),
             src->record.name);
    return false;
}

static bool length_matches(refseq_source *src, size_t file, int32_t tid)
{
    hts_pos_t expected = cm_bam_reflen(header_of(src, file), tid);

    if ((hts_pos_t)src->record.len == expected)
        return true;

    snprintf(src->error, sizeof src->error,
             "reference \"%s\" is %" PRIhts_pos " bases in the header of %s "
             "but %zu in the FASTA",
             src->record.name, expected, path_of(src, file), src->record.len);
    return false;
}

/* The MD5 of the record in hand, taken once however many files ask for it, and
 * not at all where none does. */
typedef struct {
    char value[CHECKSUM_LEN + 1];
    bool taken;
} digest;

static const char *digest_of(digest *md5, const cm_fasta_record *record)
{
    if (!md5->taken) {
        if (!checksum_sequence(record->seq, record->len, md5->value))
            return NULL;

        md5->taken = true;
    }

    return md5->value;
}

/* Whether the sequence is the one the alignments were made against, where the
 * header says what that was.
 *
 * A name and a length are a description of a reference rather than the
 * reference itself: another sequence answering to both passes them and is
 * scored against regardless. M5 is what settles it, being taken over the bases.
 * It is optional and frequently absent, though, so a reference declaring none
 * is taken on trust -- refusing those would refuse most files that exist. One
 * declaring something that is not an MD5 is a different matter: the header is
 * wrong, and going on would mean passing a check that was never made. */
static bool checksum_matches(refseq_source *src, size_t file, int32_t tid, digest *md5)
{
    size_t      declared_len = 0;
    const char *declared     = cm_bam_sq_checksum(&src->declarations[file], tid,
                                                  &declared_len);
    const char *computed;

    if (!declared)
        return true;

    if (declared_len != CHECKSUM_LEN) {
        snprintf(src->error, sizeof src->error,
                 "%s declares an M5 of %zu characters for reference \"%s\", "
                 "which is not an MD5 checksum",
                 path_of(src, file), declared_len, src->record.name);
        return false;
    }

    computed = digest_of(md5, &src->record);

    if (!computed) {
        snprintf(src->error, sizeof src->error,
                 "unable to compute a checksum for reference \"%s\"",
                 src->record.name);
        return false;
    }

    if (strncasecmp(declared, computed, CHECKSUM_LEN) == 0)
        return true;

    snprintf(src->error, sizeof src->error,
             "reference \"%s\" in the FASTA is not the sequence the alignments "
             "were made against: %s declares MD5 %.*s, the FASTA holds %s",
             src->record.name, path_of(src, file), (int)CHECKSUM_LEN, declared,
             computed);
    return false;
}

static bool matches_header(refseq_source *src, size_t file, int32_t tid, digest *md5)
{
    return name_matches(src, file, tid) &&
           length_matches(src, file, tid) &&
           checksum_matches(src, file, tid, md5);
}

static bool matches_headers(refseq_source *src, int32_t tid)
{
    digest md5 = { 0 };

    for (size_t file = 0; file < cm_bam_stream_count(src->bams); file++)
        if (!matches_header(src, file, tid, &md5))
            return false;

    return true;
}

/* ------------------------------------------------------------------------ */
/* Advance                                                                   */
/* ------------------------------------------------------------------------ */

static bool consume_through(refseq_source *src, int32_t tid)
{
    while (src->loaded <= tid) {
        int status = cm_fasta_next(src->fasta, &src->record);

        if (status == CM_ITER_EOF) {
            snprintf(src->error, sizeof src->error,
                     "FASTA ended after %d records, but the BAM header declares %d",
                     src->loaded, cm_bam_stream_nref(src->bams));
            return false;
        }

        if (status == CM_ITER_ERROR) {
            snprintf(src->error, sizeof src->error,
                     "%s", cm_fasta_error(src->fasta));
            return false;
        }

        src->loaded++;
    }

    return true;
}

const cm_fasta_record *refseq_advance(refseq_source *src, int32_t tid)
{
    if (tid < src->loaded - 1) {
        snprintf(src->error, sizeof src->error,
                 "reference %d requested after %d; the BAM is not in header order",
                 tid, src->loaded - 1);
        return NULL;
    }

    if (!consume_through(src, tid))
        return NULL;

    return matches_headers(src, tid) ? &src->record : NULL;
}

const char *refseq_error(const refseq_source *src)
{
    return src->error[0] ? src->error : NULL;
}
