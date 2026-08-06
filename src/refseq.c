/* refseq.c -- lockstep FASTA streaming, validated against the BAM header.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "refseq.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define REFSEQ_ERROR_MAX 512

struct refseq_source {
    cm_fasta_reader     *fasta;
    const cm_bam_reader *bam;
    cm_fasta_record      record;   /* the record for ordinal loaded - 1 */
    int32_t              loaded;   /* number of records consumed so far */
    char                 error[REFSEQ_ERROR_MAX];
};

refseq_source *refseq_open(const char *fasta_path, const cm_bam_reader *reader)
{
    refseq_source *src = calloc(1, sizeof *src);
    if (!src)
        return NULL;

    src->fasta = cm_fasta_open(fasta_path);
    if (!src->fasta) {
        free(src);
        return NULL;
    }

    src->bam = reader;
    return src;
}

void refseq_close(refseq_source *src)
{
    if (!src)
        return;

    cm_fasta_close(src->fasta);
    free(src);
}

/* ------------------------------------------------------------------------ */
/* Validation                                                                */
/* ------------------------------------------------------------------------ */

static bool matches_header(refseq_source *src, int32_t tid)
{
    const char *expected_name = cm_bam_refname(src->bam, tid);
    hts_pos_t   expected_len  = cm_bam_reflen(src->bam, tid);

    if (!expected_name || strcmp(src->record.name, expected_name) != 0) {
        snprintf(src->error, sizeof src->error,
                 "reference %d is \"%s\" in the BAM header but \"%s\" in the FASTA; "
                 "the two must be in the same order",
                 tid, expected_name ? expected_name : "(absent)", src->record.name);
        return false;
    }

    if ((hts_pos_t)src->record.len != expected_len) {
        snprintf(src->error, sizeof src->error,
                 "reference \"%s\" is %" PRIhts_pos " bases in the BAM header "
                 "but %zu in the FASTA",
                 src->record.name, expected_len, src->record.len);
        return false;
    }

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
                     src->loaded, cm_bam_nref(src->bam));
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

    return matches_header(src, tid) ? &src->record : NULL;
}

const char *refseq_error(const refseq_source *src)
{
    return src->error[0] ? src->error : NULL;
}
