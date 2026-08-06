/* cmuts-gen.c -- generates alignments and the reference they came from.
 *
 * Writes a BAM and a matching FASTA, for tests that need a known shape and for
 * benchmarks that need volume. Reads come out in coordinate order by
 * construction, so nothing has to be sorted afterwards however many there are.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/kstring.h>
#include <htslib/sam.h>

#include "cli.h"
#include "error.h"
#include "sample.h"
#include "simulate.h"
#include "version.h"

#define FASTA_LINE 60

enum { FORMAT_BAM, FORMAT_SAM };

/* Every parameter that varies is written as a spec -- a constant, a range, or
 * a list -- so the same option pins a value down for a test and spreads it for
 * a benchmark. They arrive as text and are parsed once, after the command line
 * is read. */
typedef struct {
    const char *output;
    int         format;

    size_t      references;
    const char *ref_length;
    double      covered;
    const char *reads_per_ref;

    const char *read_length;
    double      mismatch_rate;
    const char *insertions;
    const char *insertion_length;
    const char *deletions;
    const char *deletion_length;
    const char *soft_clips;
    const char *soft_clip_length;
    const char *mapq;
    const char *base_quality;
    double      reverse;
    const char *unmapped;

    size_t      seed;
} gen_args;

static const cli_choice FORMAT_CHOICES[] = {
    { "bam", FORMAT_BAM },
    { "sam", FORMAT_SAM },
    { NULL,  0          },
};

#define SPEC_DETAIL \
    "Written as a constant, a range LOW:HIGH, or a comma separated list."

static const cli_option OPTIONS[] = {
    { .group = "Output", .name = "output", .key = 'o', .type = OPT_STRING,
      .offset = offsetof(gen_args, output), .metavar = "PREFIX",
      .help = "write PREFIX.bam and PREFIX.fasta",
      .required = true },
    { .group = "Output", .name = "format", .type = OPT_ENUM,
      .offset = offsetof(gen_args, format), .metavar = "FORMAT",
      .help = "alignment format to write",
      .detail = "SAM is readable in a diff, which is worth having for a small "
                "fixture; BAM is what anything large should be.",
      .choices = FORMAT_CHOICES },

    { .group = "Layout", .name = "references", .type = OPT_SIZE,
      .offset = offsetof(gen_args, references), .metavar = "N",
      .help = "how many references to write",
      .minimum = 1, .maximum = 1L << 30 },
    { .group = "Layout", .name = "ref-length", .type = OPT_STRING,
      .offset = offsetof(gen_args, ref_length), .metavar = "SPEC",
      .help = "length of each reference", .detail = SPEC_DETAIL },
    { .group = "Layout", .name = "covered", .type = OPT_DOUBLE,
      .offset = offsetof(gen_args, covered), .metavar = "F",
      .help = "fraction of references receiving any reads",
      .detail = "Below one, some references get nothing at all, which is the "
                "shape of a sparse experiment and the case where most rows of "
                "the output are never written.",
      .minimum = 0, .maximum = 1 },
    { .group = "Layout", .name = "reads-per-ref", .type = OPT_STRING,
      .offset = offsetof(gen_args, reads_per_ref), .metavar = "SPEC",
      .help = "reads on each covered reference", .detail = SPEC_DETAIL },

    { .group = "Reads", .name = "read-length", .type = OPT_STRING,
      .offset = offsetof(gen_args, read_length), .metavar = "SPEC",
      .help = "reference span each read covers",
      .detail = "The span drawn from the reference. What the read finally "
                "stores differs from it: insertions and soft clips make it "
                "longer, deletions shorter. " SPEC_DETAIL },
    { .group = "Reads", .name = "mapq", .type = OPT_STRING,
      .offset = offsetof(gen_args, mapq), .metavar = "SPEC",
      .help = "mapping quality of each read", .detail = SPEC_DETAIL },
    { .group = "Reads", .name = "base-quality", .type = OPT_STRING,
      .offset = offsetof(gen_args, base_quality), .metavar = "SPEC",
      .help = "PHRED score of each base", .detail = SPEC_DETAIL },
    { .group = "Reads", .name = "reverse", .type = OPT_DOUBLE,
      .offset = offsetof(gen_args, reverse), .metavar = "F",
      .help = "fraction of reads on the reverse strand",
      .minimum = 0, .maximum = 1 },
    { .group = "Reads", .name = "unmapped", .type = OPT_STRING,
      .offset = offsetof(gen_args, unmapped), .metavar = "SPEC",
      .help = "reads aligning nowhere", .detail = SPEC_DETAIL },

    { .group = "Differences from the reference", .name = "mismatch-rate",
      .type = OPT_DOUBLE, .offset = offsetof(gen_args, mismatch_rate),
      .metavar = "F", .help = "per aligned base",
      .minimum = 0, .maximum = 1 },
    { .group = "Differences from the reference", .name = "insertions",
      .type = OPT_STRING, .offset = offsetof(gen_args, insertions),
      .metavar = "SPEC", .help = "insertion events per read",
      .detail = "Counted as events with their own length rather than drawn "
                "from a per-base rate, since a rate produces many short "
                "insertions and effectively never a long one. " SPEC_DETAIL },
    { .group = "Differences from the reference", .name = "insertion-length",
      .type = OPT_STRING, .offset = offsetof(gen_args, insertion_length),
      .metavar = "SPEC", .help = "bases per insertion", .detail = SPEC_DETAIL },
    { .group = "Differences from the reference", .name = "deletions",
      .type = OPT_STRING, .offset = offsetof(gen_args, deletions),
      .metavar = "SPEC", .help = "deletion events per read", .detail = SPEC_DETAIL },
    { .group = "Differences from the reference", .name = "deletion-length",
      .type = OPT_STRING, .offset = offsetof(gen_args, deletion_length),
      .metavar = "SPEC", .help = "bases per deletion", .detail = SPEC_DETAIL },
    { .group = "Differences from the reference", .name = "soft-clips",
      .type = OPT_STRING, .offset = offsetof(gen_args, soft_clips),
      .metavar = "SPEC", .help = "clipped ends per read, none through both",
      .detail = "Soft-clipped bases are stored but align nowhere, so they "
                "lengthen a read without lengthening its span. " SPEC_DETAIL },
    { .group = "Differences from the reference", .name = "soft-clip-length",
      .type = OPT_STRING, .offset = offsetof(gen_args, soft_clip_length),
      .metavar = "SPEC", .help = "bases per clipped end", .detail = SPEC_DETAIL },

    { .group = "Determinism", .name = "seed", .type = OPT_SIZE,
      .offset = offsetof(gen_args, seed), .metavar = "N",
      .help = "everything generated follows from this",
      .minimum = 0, .maximum = 1L << 62 },

    { .group = "Information", .name = "help", .key = 'h', .type = OPT_FLAG,
      .help = "show this help and exit", .action = CLI_SHOW_HELP },
    { .group = "Information", .name = "version", .key = 'V', .type = OPT_FLAG,
      .help = "show the version and exit", .action = CLI_SHOW_VERSION },
    { .group = "Information", .name = "dump-options", .type = OPT_FLAG,
      .help = "describe every argument as JSON and exit",
      .hidden = true, .action = CLI_DUMP_OPTIONS },
};

static gen_args gen_defaults(void)
{
    return (gen_args){
        .format           = FORMAT_BAM,
        .references       = 100,
        .ref_length       = "400",
        .covered          = 1.0,
        .reads_per_ref    = "20",
        .read_length      = "100:400",
        .mapq             = "0,1,10,30,42,60",
        .base_quality     = "30:40",
        .reverse          = 0.5,
        .unmapped         = "10",
        .mismatch_rate    = 0.01,
        .insertions       = "0:1",
        .insertion_length = "1:5",
        .deletions        = "0:1",
        .deletion_length  = "1:3",
        .soft_clips       = "0:2",
        .soft_clip_length = "5:30",
        .seed             = 1,
    };
}

/* ------------------------------------------------------------------------ */
/* Turning the text into a model                                             */
/* ------------------------------------------------------------------------ */

typedef struct {
    spec ref_length;
    spec reads_per_ref;
    spec unmapped;
} layout_specs;

static int parse_one(spec *out, const char *name, const char *text,
                     char *error, size_t error_len)
{
    /* The option's name goes in first and the reason is written after it, so
     * there is no second buffer whose contents might not fit once the two are
     * put together. */
    int prefix = snprintf(error, error_len, "--%s: ", name);

    if (prefix < 0 || (size_t)prefix >= error_len)
        return -1;

    if (spec_parse(out, text, error + prefix, error_len - (size_t)prefix) != 0)
        return -1;

    error[0] = '\0';
    return 0;
}

static int build_model(sim_model *model, layout_specs *layout, const gen_args *args,
                       char *error, size_t error_len)
{
    model->mismatch_rate    = args->mismatch_rate;
    model->reverse_fraction = args->reverse;

    return parse_one(&layout->ref_length, "ref-length", args->ref_length, error, error_len)
        || parse_one(&layout->reads_per_ref, "reads-per-ref", args->reads_per_ref, error, error_len)
        || parse_one(&layout->unmapped, "unmapped", args->unmapped, error, error_len)
        || parse_one(&model->length, "read-length", args->read_length, error, error_len)
        || parse_one(&model->mapq, "mapq", args->mapq, error, error_len)
        || parse_one(&model->base_quality, "base-quality", args->base_quality, error, error_len)
        || parse_one(&model->insertions, "insertions", args->insertions, error, error_len)
        || parse_one(&model->insertion_length, "insertion-length", args->insertion_length, error, error_len)
        || parse_one(&model->deletions, "deletions", args->deletions, error, error_len)
        || parse_one(&model->deletion_length, "deletion-length", args->deletion_length, error, error_len)
        || parse_one(&model->soft_clips, "soft-clips", args->soft_clips, error, error_len)
        || parse_one(&model->soft_clip_length, "soft-clip-length", args->soft_clip_length, error, error_len)
        ? -1 : 0;
}

/* ------------------------------------------------------------------------ */
/* References                                                                */
/* ------------------------------------------------------------------------ */

/* Each reference draws from its own stream, so what it looks like depends on
 * its index and the seed alone -- not on how many reads the references before
 * it happened to need. */
static void seed_for_reference(rng *r, size_t seed, size_t tid)
{
    rng_seed(r, (uint64_t)seed * 0x9e3779b97f4a7c15ULL + (uint64_t)tid);
}

static void reference_name(char *out, size_t size, size_t tid)
{
    snprintf(out, size, "ref%06zu", tid);
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

static sam_hdr_t *build_header(const gen_args *args, const layout_specs *layout,
                               const size_t *lengths)
{
    kstring_t  text = KS_INITIALIZE;
    sam_hdr_t *hdr;
    char       name[64];

    (void)layout;
    ksprintf(&text, "@HD\tVN:1.6\tSO:coordinate\n");

    for (size_t tid = 0; tid < args->references; tid++) {
        reference_name(name, sizeof name, tid);
        ksprintf(&text, "@SQ\tSN:%s\tLN:%zu\n", name, lengths[tid]);
    }

    hdr = sam_hdr_parse(text.l, text.s);
    free(text.s);

    return hdr;
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

typedef struct {
    samFile       *out;
    sam_hdr_t     *hdr;
    bam1_t        *rec;
    sim_scratch   *scratch;
    sim_placement *places;
    size_t         places_capacity;
    size_t         next_read;
} writer;

static int write_reference_reads(writer *w, const sim_model *model, rng *r,
                                 int32_t tid, const char *seq, size_t reflen,
                                 size_t reads)
{
    char name[64];

    if (reads > w->places_capacity)
        reads = w->places_capacity;

    for (size_t i = 0; i < reads; i++)
        w->places[i] = sim_place(model, r, reflen);

    qsort(w->places, reads, sizeof *w->places, by_start);

    for (size_t i = 0; i < reads; i++) {
        snprintf(name, sizeof name, "read%zu", w->next_read++);

        if (sim_alignment(w->rec, w->scratch, model, r, name, tid, seq, w->places[i]) < 0)
            return -1;

        if (sam_write1(w->out, w->hdr, w->rec) < 0)
            return -1;
    }

    return 0;
}

static int write_unmapped(writer *w, const sim_model *model, rng *r, size_t count)
{
    char name[64];

    for (size_t i = 0; i < count; i++) {
        snprintf(name, sizeof name, "unmapped%zu", i);

        if (sim_unmapped(w->rec, w->scratch, model, r, name) < 0 ||
            sam_write1(w->out, w->hdr, w->rec) < 0)
            return -1;
    }

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Generating                                                                */
/* ------------------------------------------------------------------------ */

static size_t *draw_reference_lengths(const gen_args *args, const layout_specs *layout)
{
    size_t *lengths = calloc(args->references, sizeof *lengths);
    rng     r;

    if (!lengths)
        return NULL;

    for (size_t tid = 0; tid < args->references; tid++) {
        long len;

        seed_for_reference(&r, args->seed, tid);
        len = spec_draw(&layout->ref_length, &r);
        lengths[tid] = len < 1 ? 1 : (size_t)len;
    }

    return lengths;
}

static int generate(const gen_args *args, const sim_model *model,
                    const layout_specs *layout, char *error, size_t error_len)
{
    char     path[4096];
    size_t  *lengths = draw_reference_lengths(args, layout);
    writer   w       = { 0 };
    FILE    *fasta   = NULL;
    char    *seq     = NULL;
    size_t   longest = 0;
    rng      r;
    int      status  = -1;

    if (!lengths) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    for (size_t tid = 0; tid < args->references; tid++)
        if (lengths[tid] > longest)
            longest = lengths[tid];

    w.hdr             = build_header(args, layout, lengths);
    w.rec             = bam_init1();
    w.scratch         = sim_scratch_create(model);
    w.places_capacity = (size_t)spec_maximum(&layout->reads_per_ref) + 1;
    w.places          = calloc(w.places_capacity, sizeof *w.places);
    seq               = malloc(longest + 1);

    snprintf(path, sizeof path, "%s.%s", args->output,
             args->format == FORMAT_SAM ? "sam" : "bam");
    w.out = sam_open(path, args->format == FORMAT_SAM ? "w" : "wb");

    snprintf(path, sizeof path, "%s.fasta", args->output);
    fasta = fopen(path, "w");

    if (!w.hdr || !w.rec || !w.scratch || !w.places || !seq || !w.out || !fasta) {
        snprintf(error, error_len, "unable to prepare the output");
        goto done;
    }

    if (sam_hdr_write(w.out, w.hdr) < 0) {
        snprintf(error, error_len, "unable to write the alignment header");
        goto done;
    }

    for (size_t tid = 0; tid < args->references; tid++) {
        char   name[64];
        size_t reflen = lengths[tid];
        size_t reads;

        seed_for_reference(&r, args->seed, tid);
        (void)spec_draw(&layout->ref_length, &r);   /* the draw already spent */

        reference_name(name, sizeof name, tid);
        sim_sequence(seq, reflen, &r);
        write_fasta_record(fasta, name, seq, reflen);

        if (!rng_chance(&r, args->covered))
            continue;

        reads = (size_t)spec_draw(&layout->reads_per_ref, &r);
        if (write_reference_reads(&w, model, &r, (int32_t)tid, seq, reflen, reads) < 0) {
            snprintf(error, error_len, "unable to write an alignment");
            goto done;
        }
    }

    rng_seed(&r, (uint64_t)args->seed ^ 0xa5a5a5a5a5a5a5a5ULL);
    if (write_unmapped(&w, model, &r, (size_t)spec_draw(&layout->unmapped, &r)) < 0) {
        snprintf(error, error_len, "unable to write an unmapped read");
        goto done;
    }

    status = 0;

done:
    if (fasta)
        fclose(fasta);
    if (w.out)
        sam_close(w.out);
    free(seq);
    free(w.places);
    sim_scratch_destroy(w.scratch);
    bam_destroy1(w.rec);
    sam_hdr_destroy(w.hdr);
    free(lengths);

    return status;
}

/* ------------------------------------------------------------------------ */
/* Entry point                                                               */
/* ------------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    gen_args     defaults = gen_defaults();
    gen_args     args;
    sim_model    model    = { 0 };
    layout_specs layout   = { 0 };
    char         error[CM_ERROR_MAX];

    cli_spec spec = {
        .program   = "cmuts-gen",
        .version   = CMUTS_VERSION,
        .summary   = "generate alignments and the reference they came from",
        .options   = OPTIONS,
        .n_options = sizeof OPTIONS / sizeof *OPTIONS,
        .defaults  = &defaults,
        .args_size = sizeof defaults,
    };

    switch (cli_parse(&spec, argc, argv, &args)) {
        case CLI_DONE:  return 0;
        case CLI_ERROR: return 2;
        case CLI_OK:    break;
    }

    if (build_model(&model, &layout, &args, error, sizeof error) < 0 ||
        generate(&args, &model, &layout, error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", spec.program, error);
        return 1;
    }

    return 0;
}
