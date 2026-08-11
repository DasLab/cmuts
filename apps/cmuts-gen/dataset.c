/* dataset.c -- writing the alignments and the reference they came from.
 *
 * A reference's reads are all placed before any of them is built, so the records can be
 * emitted in coordinate order however many there are and nothing has to be sorted afterwards.
 *
 * A reference's content derives from the seed and its index alone, so any one of them is
 * reproducible without generating those before it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "dataset.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <htslib/kstring.h>
#include <htslib/sam.h>

/* Bases per line of a FASTA record. */
#define FASTA_LINE 60

/* Names are generated, so one bound serves them all. */
#define NAME_LEN 64

#define PATH_LEN 4096

/* ------------------------------------------------------------------------ */
/* Names and paths                                                           */
/* ------------------------------------------------------------------------ */

/* Zero-padded, so that the references sort in order wherever they are listed. */
static void reference_name(char *out, size_t size, size_t tid)
{
    snprintf(out, size, "ref%06zu", tid);
}

static void read_name(char *out, size_t size, const char *kind, size_t index)
{
    snprintf(out, size, "%s%zu", kind, index);
}

/* Builds an output path, returning -1 where the prefix leaves no room for the extension rather
 * than writing to the truncated path. */
static int output_path(char *out, size_t size, const char *prefix, const char *extension)
{
    int written = snprintf(out, size, "%s.%s", prefix, extension);

    return written < 0 || (size_t)written >= size ? -1 : 0;
}

/* The format decides both the extension and the mode htslib opens with, so the two are named
 * together and neither is restated where a file is opened. */
typedef struct {
    const char *extension;
    const char *mode;
} format_encoding;

static const format_encoding ENCODINGS[] = {
    [DATASET_BAM] = { "bam", "wb" },
    [DATASET_SAM] = { "sam", "w"  },
};

/* ------------------------------------------------------------------------ */
/* Streams                                                                   */
/* ------------------------------------------------------------------------ */

/* What each reference draws from.
 *
 * One stream per reference and per purpose, so that a reference depends on the seed and its
 * index alone: not on how many reads the ones before it needed, and not on how many values
 * another purpose drew first. Sharing a stream between two purposes would make adding a draw
 * to one silently change what the other produced. */
typedef enum {
    STREAM_LENGTH,    /* how long a reference is */
    STREAM_CONTENT,   /* its sequence, and the reads laid on it */
    STREAM_UNMAPPED,  /* the reads that align nowhere */
} stream;

static void seed_stream(rng *r, size_t seed, size_t tid, stream purpose)
{
    rng_seed(r, ((uint64_t)seed * 0x9e3779b97f4a7c15ULL)
              ^ ((uint64_t)tid * 0xd1342543de82ef95ULL)
              ^ ((uint64_t)purpose * 0xa24baed4963ee407ULL));
}

/* A count drawn from a distribution, held to limit. A spec may be written with a negative
 * bound, and nothing is drawn fewer than zero times. */
static size_t draw_count(const distribution *d, rng *r, size_t limit)
{
    long drawn = distribution_draw(d, r);

    if (drawn < 0) {
        return 0;
    }

    return (size_t)drawn < limit ? (size_t)drawn : limit;
}

/* ------------------------------------------------------------------------ */
/* References                                                                */
/* ------------------------------------------------------------------------ */

/* Every reference's length, drawn up front, the header having to declare them all before the
 * first record is written. */
static size_t *draw_reference_lengths(const dataset_config *cfg)
{
    size_t *lengths = calloc(cfg->references, sizeof *lengths);

    if (!lengths) {
        return NULL;
    }

    for (size_t tid = 0; tid < cfg->references; tid++) {
        rng  r;
        long len;

        seed_stream(&r, cfg->seed, tid, STREAM_LENGTH);
        len = distribution_draw(&cfg->ref_length, &r);

        /* No read could be placed on a reference of no bases. */
        lengths[tid] = len < 1 ? 1 : (size_t)len;
    }

    return lengths;
}

static size_t longest_reference(const size_t *lengths, size_t n)
{
    size_t longest = 0;

    for (size_t i = 0; i < n; i++) {
        if (lengths[i] > longest) {
            longest = lengths[i];
        }
    }

    return longest;
}

static sam_hdr_t *build_header(const dataset_config *cfg, const size_t *lengths)
{
    kstring_t  text = KS_INITIALIZE;
    sam_hdr_t *hdr;

    ksprintf(&text, "@HD\tVN:1.6\tSO:coordinate\n");

    for (size_t tid = 0; tid < cfg->references; tid++) {
        char name[NAME_LEN];

        reference_name(name, sizeof name, tid);
        ksprintf(&text, "@SQ\tSN:%s\tLN:%zu\n", name, lengths[tid]);
    }

    hdr = sam_hdr_parse(text.l, text.s);
    free(text.s);

    return hdr;
}

static void write_fasta_record(FILE *fasta, const char *name, const char *seq, size_t len)
{
    fprintf(fasta, ">%s\n", name);

    for (size_t i = 0; i < len; i += FASTA_LINE) {
        size_t run = len - i < FASTA_LINE ? len - i : FASTA_LINE;

        fwrite(seq + i, 1, run, fasta);
        fputc('\n', fasta);
    }
}

/* ------------------------------------------------------------------------ */
/* The writer                                                                */
/* ------------------------------------------------------------------------ */

/* Both files, and everything reused from one read to the next. The buffers are sized to the
 * largest the configuration can produce, so nothing is allocated once writing has begun. */
typedef struct {
    const dataset_config *cfg;

    samFile   *alignments;
    sam_hdr_t *header;
    FILE      *fasta;

    bam1_t        *rec;
    sim_scratch   *scratch;
    sim_placement *places;    /* one reference's reads, ordered by start */
    size_t         capacity;  /* reads a reference may receive */
    char          *sequence;  /* the reference being written */

    size_t next_read;         /* numbering, over the run as a whole */
} writer;

/* Room for the most reads a reference can receive, and never for none, an empty allocation
 * being indistinguishable from a failed one. */
static size_t places_needed(const distribution *reads_per_ref)
{
    long most = distribution_maximum(reads_per_ref);

    return (most > 0 ? (size_t)most : 0) + 1;
}

static int writer_open(writer *w, const dataset_config *cfg, const size_t *lengths)
{
    const format_encoding *encoding = &ENCODINGS[cfg->format];
    char                   path[PATH_LEN];

    w->cfg      = cfg;
    w->capacity = places_needed(&cfg->reads_per_ref);

    w->header   = build_header(cfg, lengths);
    w->rec      = bam_init1();
    w->scratch  = sim_scratch_create(&cfg->model);
    w->places   = calloc(w->capacity, sizeof *w->places);
    w->sequence = malloc(longest_reference(lengths, cfg->references) + 1);

    if (output_path(path, sizeof path, cfg->prefix, encoding->extension) < 0) {
        return -1;
    }
    w->alignments = sam_open(path, encoding->mode);

    if (output_path(path, sizeof path, cfg->prefix, "fasta") < 0) {
        return -1;
    }
    w->fasta = fopen(path, "w");

    return w->header && w->rec && w->scratch && w->places && w->sequence &&
           w->alignments && w->fasta ? 0 : -1;
}

static void writer_close(writer *w)
{
    if (w->fasta) {
        fclose(w->fasta);
    }
    if (w->alignments) {
        sam_close(w->alignments);
    }

    free(w->sequence);
    free(w->places);
    sim_scratch_destroy(w->scratch);
    bam_destroy1(w->rec);
    sam_hdr_destroy(w->header);
}

/* ------------------------------------------------------------------------ */
/* Reads                                                                     */
/* ------------------------------------------------------------------------ */

static int by_start(const void *a, const void *b)
{
    hts_pos_t x = ((const sim_placement *)a)->start;
    hts_pos_t y = ((const sim_placement *)b)->start;

    return x < y ? -1 : x > y;
}

static int write_read(writer *w, rng *r, int32_t tid, sim_placement where)
{
    char name[NAME_LEN];

    read_name(name, sizeof name, "read", w->next_read++);

    if (sim_alignment(w->rec, w->scratch, &w->cfg->model, r, name, tid,
                      w->sequence, where) < 0) {
        return -1;
    }

    return sam_write1(w->alignments, w->header, w->rec) < 0 ? -1 : 0;
}

/* Writes every read of one reference. The placements are eight bytes apiece and are drawn and
 * sorted before any record is built, so coordinate order costs the sort rather than holding
 * every record of the reference in memory. */
static int write_reference_reads(writer *w, rng *r, int32_t tid, size_t reflen)
{
    size_t reads = draw_count(&w->cfg->reads_per_ref, r, w->capacity);

    for (size_t i = 0; i < reads; i++) {
        w->places[i] = sim_place(&w->cfg->model, r, reflen);
    }

    qsort(w->places, reads, sizeof *w->places, by_start);

    for (size_t i = 0; i < reads; i++) {
        if (write_read(w, r, tid, w->places[i]) < 0) {
            return -1;
        }
    }

    return 0;
}

/* Writes one reference and its reads. The sequence is written whether or not the reference
 * receives any reads, a reference with none being a row the output holds and nothing in it --
 * which is what a sparse experiment consists of. */
static int write_reference(writer *w, size_t tid, size_t reflen)
{
    char name[NAME_LEN];
    rng  r;

    seed_stream(&r, w->cfg->seed, tid, STREAM_CONTENT);
    reference_name(name, sizeof name, tid);

    sim_sequence(w->sequence, reflen, &r);
    write_fasta_record(w->fasta, name, w->sequence, reflen);

    if (!rng_chance(&r, w->cfg->covered)) {
        return 0;
    }

    return write_reference_reads(w, &r, (int32_t)tid, reflen);
}

/* Writes the reads belonging to no reference. Drawn once for the run and written after
 * everything that is placed, where a coordinate-sorted file keeps them. */
static int write_unmapped_reads(writer *w)
{
    rng    r;
    size_t count;

    seed_stream(&r, w->cfg->seed, 0, STREAM_UNMAPPED);
    count = draw_count(&w->cfg->unmapped, &r, SIZE_MAX);

    for (size_t i = 0; i < count; i++) {
        char name[NAME_LEN];

        read_name(name, sizeof name, "unmapped", i);

        if (sim_unmapped(w->rec, w->scratch, &w->cfg->model, &r, name) < 0 ||
            sam_write1(w->alignments, w->header, w->rec) < 0) {
            return -1;
        }
    }

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Assembly                                                                  */
/* ------------------------------------------------------------------------ */

int dataset_write(const dataset_config *cfg, char *error, size_t error_len)
{
    writer  w       = { 0 };
    size_t *lengths = draw_reference_lengths(cfg);
    int     status  = -1;

    if (!lengths) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    if (writer_open(&w, cfg, lengths) < 0) {
        snprintf(error, error_len, "unable to prepare the output");
        goto done;
    }

    if (sam_hdr_write(w.alignments, w.header) < 0) {
        snprintf(error, error_len, "unable to write the alignment header");
        goto done;
    }

    for (size_t tid = 0; tid < cfg->references; tid++) {
        if (write_reference(&w, tid, lengths[tid]) < 0) {
            snprintf(error, error_len, "unable to write an alignment");
            goto done;
        }
    }

    if (write_unmapped_reads(&w) < 0) {
        snprintf(error, error_len, "unable to write an unmapped read");
        goto done;
    }

    status = 0;

done:
    writer_close(&w);
    free(lengths);

    return status;
}
