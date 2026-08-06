/* bam.c -- BAM iteration, backed by htslib's sam reader.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "bam.h"

#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <htslib/bgzf.h>
#include <htslib/cram.h>
#include <htslib/hfile.h>


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
    uint64_t    header_end;  /* where the alignments begin */
    uint64_t    file_size;   /* 0 where the input is not a file of known size */
    bool        at_end;
    const char *error;
};

/* ------------------------------------------------------------------------ */
/* Position                                                                  */
/* ------------------------------------------------------------------------ */

/* How far into the file the reader has read, in compressed bytes. Every format
 * keeps the figure somewhere else, and none of them keeps a record count. */
static uint64_t raw_offset(const cm_bam_reader *reader)
{
    const htsFile *file = reader->file;

    if (file->is_cram)
        return (uint64_t)htell(cram_fd_get_fp(file->fp.cram));

    /* A bgzipped SAM reads through BGZF as well, so the format alone does not
     * settle which of these applies. */
    if (file->is_bgzf)
        return (uint64_t)(bgzf_tell(file->fp.bgzf) >> 16);

    return file->fp.hfile ? (uint64_t)htell(file->fp.hfile) : 0;
}

static uint64_t size_of(const char *path)
{
    struct stat info;

    return stat(path, &info) == 0 && info.st_size > 0 ? (uint64_t)info.st_size : 0;
}

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

    reader->header_end = raw_offset(reader);
    reader->file_size  = size_of(path);

    return reader;
}

void cm_bam_close(cm_bam_reader *reader)
{
    if (reader)
        reader_free(reader);
}

int cm_bam_set_reference(cm_bam_reader *reader, const char *fasta_path)
{
    if (hts_set_opt(reader->file, CRAM_OPT_REFERENCE, fasta_path) < 0) {
        reader->error = "unable to point the reader at the reference";
        return -1;
    }

    return 0;
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

    if (status == SAM_END_OF_FILE) {
        reader->at_end = true;
        return CM_ITER_EOF;
    }

    reader->error = "error reading alignment record";
    return CM_ITER_ERROR;
}

const char *cm_bam_error(const cm_bam_reader *reader)
{
    return reader->error;
}

/* An offset points at the block being read rather than past it, so it stops
 * short of the end even once there is nothing left to read. */
uint64_t cm_bam_position(const cm_bam_reader *reader)
{
    uint64_t at;

    if (reader->at_end)
        return cm_bam_span(reader);

    at = raw_offset(reader);

    return at > reader->header_end ? at - reader->header_end : 0;
}

uint64_t cm_bam_span(const cm_bam_reader *reader)
{
    return reader->file_size > reader->header_end
         ? reader->file_size - reader->header_end
         : 0;
}

/* ------------------------------------------------------------------------ */
/* Header text                                                               */
/* ------------------------------------------------------------------------ */

/* The header is read as the text it arrived as. Asking htslib for a tag would
 * have it parse the whole header into records first, and on a file with
 * millions of references that costs gigabytes to answer a question about a few
 * characters -- 3.85 GB against a header text of 0.50 GB, measured on a file of
 * 24 million. The text is already in memory and holds the same information. */

static const char *line_end(const char *line)
{
    const char *newline = strchr(line, '\n');

    return newline ? newline : line + strlen(line);
}

static const char *next_line(const char *line)
{
    const char *newline = strchr(line, '\n');

    return (newline && newline[1]) ? newline + 1 : NULL;
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

/* Where the field a value belongs to ends. */
static const char *field_end(const char *value, const char *end)
{
    const char *tab = memchr(value, '\t', (size_t)(end - value));

    return tab ? tab : end;
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

/* @HD is the first line where a sort order appears at all, so only the first
 * line is examined. A header carrying one elsewhere violates the specification
 * and is treated as saying nothing about sort order. */
bool cm_bam_is_coordinate_sorted(const cm_bam_reader *reader)
{
    const char *text = sam_hdr_str(reader->header);
    const char *end;
    const char *order;
    size_t      len = strlen(SORT_ORDER_COORDINATE);

    if (!text || strncmp(text, "@HD", 3) != 0)
        return false;

    end   = line_end(text);
    order = find_tag(text, end, "SO:");

    if (!order)
        return false;

    /* The value runs to the end of its field, so a longer word beginning with
     * "coordinate" is not one. */
    return (size_t)(field_end(order, end) - order) == len &&
           strncmp(order, SORT_ORDER_COORDINATE, len) == 0;
}

/* ------------------------------------------------------------------------ */
/* Reference declarations                                                    */
/* ------------------------------------------------------------------------ */

/* The first @SQ line at or after the given one, or NULL where there is none. */
static const char *next_declaration(const char *line)
{
    for (; line; line = next_line(line))
        if (strncmp(line, "@SQ", 3) == 0)
            return line;

    return NULL;
}

void cm_bam_sq_open(cm_bam_sq_cursor *cursor, const cm_bam_reader *reader)
{
    const char *text = sam_hdr_str(reader->header);

    cursor->line = text ? next_declaration(text) : NULL;
    cursor->tid  = 0;
}

/* The lines appear in the order the references are numbered, so reaching one
 * is a matter of walking forward over those before it. The cursor is left on
 * the line it reached rather than past it, so that asking twice for the same
 * reference answers the same both times. */
const char *cm_bam_sq_checksum(cm_bam_sq_cursor *cursor, int32_t tid, size_t *len)
{
    const char *end;
    const char *checksum;

    while (cursor->line && cursor->tid < tid) {
        cursor->line = next_declaration(next_line(cursor->line));
        cursor->tid++;
    }

    if (!cursor->line || cursor->tid != tid)
        return NULL;

    end      = line_end(cursor->line);
    checksum = find_tag(cursor->line, end, "M5:");

    if (checksum)
        *len = (size_t)(field_end(checksum, end) - checksum);

    return checksum;
}
