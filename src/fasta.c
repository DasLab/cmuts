/* fasta.c -- FASTA iteration, backed by htslib's kseq parser.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "fasta.h"

#include <stdlib.h>

#include <htslib/bgzf.h>
#include <htslib/kseq.h>

KSEQ_INIT(BGZF *, bgzf_read)

/* kseq_read() reports these; kseq.h documents but does not name them. A
 * non-negative return is the length of the sequence just parsed. */
enum {
    KSEQ_END_OF_FILE       = -1,
    KSEQ_TRUNCATED_QUALITY = -2,
    KSEQ_STREAM_ERROR      = -3,
    KSEQ_OVERFLOW          = -4,
};

struct cm_fasta_reader {
    BGZF       *file;
    kseq_t     *seq;
    const char *error;
};

/* ------------------------------------------------------------------------ */
/* Lifetime                                                                  */
/* ------------------------------------------------------------------------ */

static void reader_free(cm_fasta_reader *reader)
{
    if (reader->seq)
        kseq_destroy(reader->seq);
    if (reader->file)
        bgzf_close(reader->file);
    free(reader);
}

cm_fasta_reader *cm_fasta_open(const char *path)
{
    cm_fasta_reader *reader = calloc(1, sizeof *reader);
    if (!reader)
        return NULL;

    reader->file = bgzf_open(path, "r");
    if (!reader->file) {
        reader_free(reader);
        return NULL;
    }

    reader->seq = kseq_init(reader->file);
    if (!reader->seq) {
        reader_free(reader);
        return NULL;
    }

    return reader;
}

void cm_fasta_close(cm_fasta_reader *reader)
{
    if (reader)
        reader_free(reader);
}

/* ------------------------------------------------------------------------ */
/* Iteration                                                                 */
/* ------------------------------------------------------------------------ */

/* kseq leaves the comment empty rather than absent when a header line carries
 * no description. */
static const char *record_comment(const kseq_t *seq)
{
    return seq->comment.l > 0 ? seq->comment.s : NULL;
}

static void record_from_kseq(const kseq_t *seq, cm_fasta_record *out)
{
    out->name    = seq->name.s;
    out->comment = record_comment(seq);
    out->seq     = seq->seq.s;
    out->len     = seq->seq.l;
}

static const char *parse_error_message(int status)
{
    switch (status) {
        case KSEQ_TRUNCATED_QUALITY: return "truncated quality string";
        case KSEQ_STREAM_ERROR:      return "error reading stream";
        case KSEQ_OVERFLOW:          return "record too long";
        default:                     return "malformed record";
    }
}

int cm_fasta_next(cm_fasta_reader *reader, cm_fasta_record *out)
{
    int status = kseq_read(reader->seq);

    if (status >= 0) {
        record_from_kseq(reader->seq, out);
        return CM_ITER_OK;
    }

    if (status == KSEQ_END_OF_FILE)
        return CM_ITER_EOF;

    reader->error = parse_error_message(status);
    return CM_ITER_ERROR;
}

const char *cm_fasta_error(const cm_fasta_reader *reader)
{
    return reader->error;
}
