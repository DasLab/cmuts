/* reads.c -- reading a file of sequencing reads, and writing one out as FASTQ.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "reads.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

/* The offset that a FASTQ adds to a PHRED score to make it printable. */
#define PHRED_OFFSET 33

/* The quality that htslib stores for every base of a record that has no qualities. */
#define QUALITY_ABSENT 0xFF

static int fail(const char *path, const char *why, char *error, size_t error_len)
{
    snprintf(error, error_len, "%s: %s", path, why);

    return -1;
}

/* ------------------------------------------------------------------------ */
/* What a file holds                                                         */
/* ------------------------------------------------------------------------ */

reads_format reads_format_of(const char *path)
{
    htsFile *file = hts_open(path, "r");

    if (!file) {
        return READS_NEITHER;
    }

    reads_format found = READS_NEITHER;

    switch (hts_get_format(file)->format) {
        case fastq_format: found = READS_FASTQ; break;
        case bam:          found = READS_BAM;   break;
        default:                                break;
    }

    hts_close(file);

    return found;
}

/* Reads the SAM flag of the first record. Sets the flag to zero if the file holds no
 * records. */
static int first_flag(samFile *file, sam_hdr_t *header, uint16_t *out)
{
    bam1_t *record = bam_init1();

    if (!record) {
        return -1;
    }

    int status = sam_read1(file, header, record);

    *out = status >= 0 ? record->core.flag : 0;
    bam_destroy1(record);

    return status >= -1 ? 0 : -1;
}

int reads_inspect_bam(const char *path, bool *aligned, bool *paired,
                      char *error, size_t error_len)
{
    samFile *file = sam_open(path, "r");

    if (!file) {
        return fail(path, "cannot be opened", error, error_len);
    }

    sam_hdr_t *header = sam_hdr_read(file);
    uint16_t   flag   = 0;
    int        status = 0;

    if (!header) {
        status = fail(path, "has no header", error, error_len);
    } else if (first_flag(file, header, &flag) < 0) {
        status = fail(path, "its records cannot be read", error, error_len);
    } else {
        *aligned = sam_hdr_nref(header) > 0;
        *paired  = (flag & BAM_FPAIRED) != 0;
    }

    sam_hdr_destroy(header);
    sam_close(file);

    return status;
}

/* ------------------------------------------------------------------------ */
/* Writing a BAM file out as FASTQ                                           */
/* ------------------------------------------------------------------------ */

/* Printable text for the bases and the qualities of one record. Both buffers grow to hold
 * the longest record that has been written. */
typedef struct {
    char  *bases;
    char  *qualities;
    size_t capacity;
} record_text;

static bool text_reserve(record_text *text, size_t len)
{
    if (len + 1 <= text->capacity) {
        return true;
    }

    char *bases     = realloc(text->bases, len + 1);
    char *qualities = realloc(text->qualities, len + 1);

    if (bases) {
        text->bases = bases;
    }

    if (qualities) {
        text->qualities = qualities;
    }

    if (!bases || !qualities) {
        return false;
    }

    text->capacity = len + 1;

    return true;
}

static void text_free(record_text *text)
{
    free(text->bases);
    free(text->qualities);
}

static void fill_text(record_text *text, const bam1_t *record)
{
    const uint8_t *bases     = bam_get_seq(record);
    const uint8_t *qualities = bam_get_qual(record);
    size_t         len       = (size_t)record->core.l_qseq;

    for (size_t i = 0; i < len; i++) {
        text->bases[i]     = seq_nt16_str[bam_seqi(bases, i)];
        text->qualities[i] = (char)(qualities[i] + PHRED_OFFSET);
    }

    text->bases[len]     = '\0';
    text->qualities[len] = '\0';
}

/* Returns why the record cannot be written as FASTQ, or NULL if it can be. An unaligned
 * BAM file holds complete reads only. */
static const char *refusal_reason(const bam1_t *record)
{
    if (record->core.l_qseq == 0) {
        return "has no sequence";
    }

    if (bam_get_qual(record)[0] == QUALITY_ABSENT) {
        return "has no base qualities";
    }

    if (record->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
        return "is a secondary or a supplementary record";
    }

    /* An unaligned read lies on no strand. A reverse flag would mean the bases are
     * stored the other way round from how the instrument called them. */
    if (record->core.flag & BAM_FREVERSE) {
        return "is marked as reverse, which an unaligned read cannot be";
    }

    return NULL;
}

static int write_record(FILE *out, record_text *text, const bam1_t *record,
                        char *error, size_t error_len)
{
    const char *why = refusal_reason(record);

    if (why) {
        snprintf(error, error_len, "read \"%s\" %s", bam_get_qname(record), why);
        return -1;
    }

    if (!text_reserve(text, (size_t)record->core.l_qseq)) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    fill_text(text, record);
    fprintf(out, "@%s\n%s\n+\n%s\n", bam_get_qname(record), text->bases, text->qualities);

    return 0;
}

/* Writes every record of an open file. The caller must read the header first. */
static int write_records(samFile *file, sam_hdr_t *header, FILE *out, const char *path,
                         char *error, size_t error_len)
{
    record_text text    = { 0 };
    bam1_t     *record  = bam_init1();
    int         status  = 0;
    int         outcome = 0;

    if (!record) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    while (status == 0 && (outcome = sam_read1(file, header, record)) >= 0) {
        status = write_record(out, &text, record, error, error_len);
    }

    /* sam_read1 returns -1 at the end of the file, and less than that where a record
     * could not be read. Stopping early on the second would align the reads that were
     * read and report nothing about the ones that were not. */
    if (status == 0 && outcome < -1) {
        status = fail(path, "could not be read to the end", error, error_len);
    }

    bam_destroy1(record);
    text_free(&text);

    return status;
}

int reads_write_fastq(const char *path, int fd, char *error, size_t error_len)
{
    FILE *out = fdopen(fd, "w");

    if (!out) {
        close(fd);
        return fail(path, "could not be passed to the aligner", error, error_len);
    }

    samFile *file = sam_open(path, "r");

    if (!file) {
        fclose(out);
        return fail(path, "cannot be opened", error, error_len);
    }

    sam_hdr_t *header = sam_hdr_read(file);
    int        status = header ? write_records(file, header, out, path, error, error_len)
                               : fail(path, "has no header", error, error_len);

    /* The descriptor is closed even after a failure, because the aligner reads until the
     * stream ends. */
    if (fclose(out) != 0 && status == 0) {
        status = fail(path, "could not be passed to the aligner", error, error_len);
    }

    sam_hdr_destroy(header);
    sam_close(file);

    return status;
}
