/* bam.c -- BAM iteration, backed by htslib's sam reader.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "bam.h"

#include <stdlib.h>
#include <string.h>


/* sam_read1() returns a non-negative value per record and this on a clean end
 * of stream; anything lower is a read error. */
#define SAM_END_OF_FILE (-1)

static const char TAG_MD[2] = { 'M', 'D' };

/* SAM spells "no qualities" as a QUAL field of "*", which BAM stores as a
 * leading 0xff in an otherwise full-length quality array. */
#define QUAL_ABSENT 0xff

static const char SORT_ORDER_COORDINATE[] = "coordinate";

struct cm_bam_reader {
    samFile    *file;
    sam_hdr_t  *header;
    bam1_t     *record;
    const char *error;
};

/* ------------------------------------------------------------------------ */
/* Lifetime                                                                  */
/* ------------------------------------------------------------------------ */

static void reader_free(cm_bam_reader *reader)
{
    if (reader->record)
        bam_destroy1(reader->record);
    if (reader->header)
        sam_hdr_destroy(reader->header);
    if (reader->file)
        sam_close(reader->file);
    free(reader);
}

cm_bam_reader *cm_bam_open(const char *path)
{
    cm_bam_reader *reader = calloc(1, sizeof *reader);
    if (!reader)
        return NULL;

    reader->file = sam_open(path, "r");
    if (!reader->file) {
        reader_free(reader);
        return NULL;
    }

    reader->header = sam_hdr_read(reader->file);
    if (!reader->header) {
        reader_free(reader);
        return NULL;
    }

    reader->record = bam_init1();
    if (!reader->record) {
        reader_free(reader);
        return NULL;
    }

    return reader;
}

void cm_bam_close(cm_bam_reader *reader)
{
    if (reader)
        reader_free(reader);
}

int cm_bam_set_threads(cm_bam_reader *reader, int threads)
{
    if (threads < 1)
        return 0;

    if (hts_set_threads(reader->file, threads) < 0) {
        reader->error = "unable to start decompression threads";
        return -1;
    }

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Field extraction                                                          */
/* ------------------------------------------------------------------------ */

static const char *md_tag(const bam1_t *record)
{
    const uint8_t *tag = bam_aux_get(record, TAG_MD);
    return tag ? bam_aux2Z(tag) : NULL;
}

static const uint8_t *quality_scores(const bam1_t *record)
{
    const uint8_t *qual = bam_get_qual(record);

    if (record->core.l_qseq <= 0 || qual[0] == QUAL_ABSENT)
        return NULL;

    return qual;
}

void cm_bam_record_view(const bam1_t *record, cm_bam_record *out)
{
    out->tid     = record->core.tid;
    out->pos     = record->core.pos;
    out->flag    = record->core.flag;
    out->mapq    = record->core.qual;
    out->cigar   = bam_get_cigar(record);
    out->n_cigar = record->core.n_cigar;
    out->md      = md_tag(record);
    out->qual    = quality_scores(record);
    out->l_qseq  = record->core.l_qseq;
}

const bam1_t *cm_bam_raw(const cm_bam_reader *reader)
{
    return reader->record;
}

/* ------------------------------------------------------------------------ */
/* Iteration                                                                 */
/* ------------------------------------------------------------------------ */

int cm_bam_next(cm_bam_reader *reader, cm_bam_record *out)
{
    int status = sam_read1(reader->file, reader->header, reader->record);

    if (status >= 0) {
        cm_bam_record_view(reader->record, out);
        return CM_ITER_OK;
    }

    if (status == SAM_END_OF_FILE)
        return CM_ITER_EOF;

    reader->error = "error reading alignment record";
    return CM_ITER_ERROR;
}

const char *cm_bam_error(const cm_bam_reader *reader)
{
    return reader->error;
}

/* ------------------------------------------------------------------------ */
/* Header queries                                                            */
/* ------------------------------------------------------------------------ */

const char *cm_bam_refname(const cm_bam_reader *reader, int32_t tid)
{
    return sam_hdr_tid2name(reader->header, tid);
}

hts_pos_t cm_bam_reflen(const cm_bam_reader *reader, int32_t tid)
{
    return sam_hdr_tid2len(reader->header, tid);
}

hts_pos_t cm_bam_max_reflen(const cm_bam_reader *reader)
{
    hts_pos_t longest = 0;
    int32_t   n       = sam_hdr_nref(reader->header);

    for (int32_t tid = 0; tid < n; tid++) {
        hts_pos_t len = sam_hdr_tid2len(reader->header, tid);
        if (len > longest)
            longest = len;
    }

    return longest;
}

int32_t cm_bam_nref(const cm_bam_reader *reader)
{
    return sam_hdr_nref(reader->header);
}

/* Where a tab-separated field with the given prefix begins, within one line. */
static const char *find_tag(const char *line, const char *end, const char *tag)
{
    size_t len = strlen(tag);

    for (const char *p = line; p + len <= end; p++)
        if (*p == '\t' && strncmp(p + 1, tag, len) == 0)
            return p + 1 + len;

    return NULL;
}

/* Asking htslib for the tag would have it parse the whole header into records
 * first, and on a file with millions of references that costs gigabytes to
 * answer a question about ten characters -- 3.85 GB against a header text of
 * 0.50 GB, measured on a file of 24 million. The text is already in memory and
 * holds the same information, so it is read directly.
 *
 * @HD is the first line where it appears at all, so only the first line is
 * examined. A header carrying it elsewhere violates the specification and is
 * treated as not saying anything about sort order. */
bool cm_bam_is_coordinate_sorted(const cm_bam_reader *reader)
{
    const char *text = sam_hdr_str(reader->header);
    const char *end;
    const char *order;
    size_t      len = strlen(SORT_ORDER_COORDINATE);

    if (!text || strncmp(text, "@HD", 3) != 0)
        return false;

    end   = strchr(text, '\n');
    end   = end ? end : text + strlen(text);
    order = find_tag(text, end, "SO:");

    if (!order || order + len > end)
        return false;

    /* The value runs to the end of the field, so a longer word starting with
     * "coordinate" is not one. */
    return strncmp(order, SORT_ORDER_COORDINATE, len) == 0 &&
           (order + len == end || order[len] == '\t');
}
