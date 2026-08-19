/* score.c -- reactivity measured against a known structure.
 *
 * A reagent reports on the bases a structure leaves open, so a reactivity should rank an
 * unpaired base above a paired one. Every metric here measures how well it does that,
 * within one reference: reactivity carries the depth and the scale of the sample it came
 * from, so values from two references do not compare.
 *
 * The HDF5 file holds no names. Its rows follow the FASTA that produced it, so the FASTA
 * names them and gives the sequence each is scored on.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "score.h"

#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "fasta.h"
#include "h5reader.h"
#include "output.h"

/* What a dot bracket record may hold. A base is paired under any bracket, open under a
 * dot, and unknown under anything else, which is left out of the scoring. */
#define BRACKETS "()[]{}<>"
#define UNPAIRED '.'

/* ------------------------------------------------------------------------ */
/* The structures                                                            */
/* ------------------------------------------------------------------------ */

/* One reference's structure. Every string owns its memory. The sequence is held only
 * where the record carried one, and is checked against the reference. */
typedef struct {
    char  *name;
    char  *seq;
    char  *pairing;
    size_t len;
} structure;

/* Every structure read, sorted by name so that a reference is found by search. */
typedef struct {
    structure *at;
    size_t     n;
    size_t     cap;
} structures;

static void structures_free(structures *set)
{
    for (size_t i = 0; i < set->n; i++) {
        free(set->at[i].name);
        free(set->at[i].seq);
        free(set->at[i].pairing);
    }

    free(set->at);
}

static int by_name(const void *a, const void *b)
{
    return strcmp(((const structure *)a)->name, ((const structure *)b)->name);
}

static const structure *structure_named(const structures *set, const char *name)
{
    /* bsearch takes a key of the element's type, so the name is copied into one. The
     * copy is read from and never written through. */
    structure key;

    memcpy(&key.name, &name, sizeof key.name);

    return bsearch(&key, set->at, set->n, sizeof *set->at, by_name);
}

static int structures_grow(structures *set)
{
    size_t     cap = set->cap ? set->cap * 2 : 64;
    structure *at  = realloc(set->at, cap * sizeof *at);

    if (!at) {
        return -1;
    }

    set->at  = at;
    set->cap = cap;

    return 0;
}

/* Returns whether the line is a pairing and not a sequence. A sequence never carries a
 * bracket or a dot, so one character decides it. */
static bool is_pairing(const char *line)
{
    return line[strspn(line, BRACKETS ".-")] == '\0' && *line != '\0';
}

/* Cuts the line at its first line ending. */
static void trim(char *line)
{
    line[strcspn(line, "\r\n")] = '\0';
}

/* Keeps one record, copying it, since the lines it came from are reused. The sequence is
 * optional and NULL where the record carried none. */
static int structures_add(structures *set, const char *name, const char *seq,
                          const char *pairing)
{
    structure *at;

    if (set->n == set->cap && structures_grow(set) < 0) {
        return -1;
    }

    at          = &set->at[set->n];
    at->name    = strdup(name);
    at->seq     = seq ? strdup(seq) : NULL;
    at->pairing = strdup(pairing);
    at->len     = strlen(pairing);

    if (!at->name || !at->pairing || (seq && !at->seq)) {
        free(at->name);
        free(at->seq);
        free(at->pairing);
        return -1;
    }

    set->n++;

    return 0;
}

/* Reads a file of dot bracket records: a header naming the reference, then the pairing,
 * with a sequence line before it where the file carries one. */
static int structures_read(structures *set, const char *path, char *error,
                           size_t error_len)
{
    char  line[4096];
    char  name[256] = { 0 };
    char *seq       = NULL;
    FILE *file      = fopen(path, "r");
    int   number    = 0;

    if (!file) {
        snprintf(error, error_len, "%s: cannot be opened", path);
        return -1;
    }

    while (fgets(line, sizeof line, file)) {
        number++;
        trim(line);

        if (line[0] == '>') {
            snprintf(name, sizeof name, "%s", line + 1 + strspn(line + 1, " \t"));
            name[strcspn(name, " \t")] = '\0';
            free(seq);
            seq = NULL;
            continue;
        }

        if (line[0] == '\0') {
            continue;
        }

        /* a record may carry its sequence before the pairing, which is kept so that the
         * reference can be checked against it */
        if (!is_pairing(line)) {
            free(seq);
            seq = strdup(line);

            if (!seq) {
                snprintf(error, error_len, "%s: out of memory", path);
                fclose(file);
                return -1;
            }

            continue;
        }

        if (name[0] == '\0') {
            snprintf(error, error_len, "%s:%d: a pairing before any name", path, number);
            free(seq);
            fclose(file);
            return -1;
        }

        if (structures_add(set, name, seq, line) < 0) {
            snprintf(error, error_len, "%s: out of memory", path);
            free(seq);
            fclose(file);
            return -1;
        }

        name[0] = '\0';
        free(seq);
        seq = NULL;
    }

    free(seq);
    fclose(file);
    qsort(set->at, set->n, sizeof *set->at, by_name);

    return 0;
}

/* ------------------------------------------------------------------------ */
/* The metrics                                                               */
/* ------------------------------------------------------------------------ */

/* One position kept for scoring. */
typedef struct {
    double value;
    bool   open;   /* whether the structure leaves the base unpaired */
} point;

/* What one reference scored. */
typedef struct {
    size_t paired;
    size_t unpaired;
    double auroc;
    double auprc;
    double mean_paired;
    double mean_unpaired;
} result;

static int by_value(const void *a, const void *b)
{
    double left  = ((const point *)a)->value;
    double right = ((const point *)b)->value;

    return left < right ? -1 : left > right ? 1 : 0;
}

/* Returns the chance an unpaired base outranks a paired one, with ties counted as half.
 * The points must be sorted by value. */
static double auroc_of(const point *at, size_t n, size_t unpaired)
{
    double sum = 0.0;
    size_t i   = 0;

    while (i < n) {
        size_t j = i;
        double rank;
        size_t open = 0;

        while (j < n && at[j].value == at[i].value) {
            open += at[j].open ? 1 : 0;
            j++;
        }

        /* every tied value takes the middle of the ranks the run covers */
        rank = ((double)i + 1.0 + (double)j) / 2.0;
        sum += rank * (double)open;
        i = j;
    }

    sum -= (double)unpaired * ((double)unpaired + 1.0) / 2.0;

    return sum / ((double)unpaired * (double)(n - unpaired));
}

/* Returns the average precision, taking the unpaired bases as what is sought. The points
 * must be sorted by value, so the walk runs from the highest reactivity down. */
static double auprc_of(const point *at, size_t n, size_t unpaired)
{
    double found = 0.0;
    double seen  = 0.0;
    double sum   = 0.0;

    for (size_t i = n; i-- > 0; ) {
        seen += 1.0;

        if (at[i].open) {
            found += 1.0;
            sum   += found / seen;
        }
    }

    return sum / (double)unpaired;
}

static double mean_of(const point *at, size_t n, bool open)
{
    double sum   = 0.0;
    size_t count = 0;

    for (size_t i = 0; i < n; i++) {
        if (at[i].open == open) {
            sum += at[i].value;
            count++;
        }
    }

    return count ? sum / (double)count : (double)NAN;
}

/* Measures the points kept for one reference. Valid only where both classes are held. */
static result measure(point *at, size_t n, size_t unpaired)
{
    qsort(at, n, sizeof *at, by_value);

    return (result){
        .paired        = n - unpaired,
        .unpaired      = unpaired,
        .auroc         = auroc_of(at, n, unpaired),
        .auprc         = auprc_of(at, n, unpaired),
        .mean_paired   = mean_of(at, n, false),
        .mean_unpaired = mean_of(at, n, true),
    };
}

/* ------------------------------------------------------------------------ */
/* One reference                                                             */
/* ------------------------------------------------------------------------ */

/* Everything one run works from. */
typedef struct {
    const score_config *cfg;
    h5reader           *reader;
    float              *reactivity;
    float              *coverage;
    point              *points;
    size_t              cap;
} context;

/* Returns whether two bases are the same, reading U and T as one another and ignoring
 * case, so a structure written as RNA matches a reference written as DNA. */
static bool same_base(char a, char b)
{
    a = (char)toupper((unsigned char)a);
    b = (char)toupper((unsigned char)b);

    a = a == 'U' ? 'T' : a;
    b = b == 'U' ? 'T' : b;

    return a == b;
}

/* Returns the first position where the record's sequence differs from the reference, or
 * -1 where they agree over the shorter of the two. */
static long disagreement(const char *seq, const cm_fasta_record *ref)
{
    for (size_t i = 0; i < ref->len && seq[i] != '\0'; i++) {
        if (!same_base(seq[i], ref->seq[i])) {
            return (long)i;
        }
    }

    return -1;
}

/* Returns the bit standing for a base, or 0 for a base the alphabet does not name, such
 * as an ambiguous one. T and U give the same bit. */
static int bit_of(char base)
{
    switch (toupper((unsigned char)base)) {
        case 'A': return SCORE_BASE_A;
        case 'C': return SCORE_BASE_C;
        case 'G': return SCORE_BASE_G;
        case 'T':
        case 'U': return SCORE_BASE_U;
        default:  return 0;
    }
}

/* Returns whether the base at this position is one the reagent modifies. */
static bool base_wanted(const context *ctx, char base)
{
    return (ctx->cfg->bases & bit_of(base)) != 0;
}

/* Fills the points from one reference's row, and gives how many are unpaired. Returns the
 * count kept. */
static size_t collect(context *ctx, const cm_fasta_record *ref,
                      const structure *known, size_t *unpaired)
{
    size_t n = 0;

    *unpaired = 0;

    for (size_t i = 0; i < ref->len && i < known->len; i++) {
        char   mark  = known->pairing[i];
        bool   open  = mark == UNPAIRED;
        double value = (double)ctx->reactivity[i];

        if (!open && !strchr(BRACKETS, mark)) {
            continue;
        }
        if (!base_wanted(ctx, ref->seq[i]) || !isfinite(value)) {
            continue;
        }
        if ((double)ctx->coverage[i] < ctx->cfg->min_coverage) {
            continue;
        }

        ctx->points[n].value = value;
        ctx->points[n].open  = open;
        *unpaired += open ? 1 : 0;
        n++;
    }

    return n;
}

/* Reads one reference's row into the context. */
static int read_row(context *ctx, int32_t tid, char *error, size_t error_len)
{
    if (h5reader_field(ctx->reader, OUT_REACTIVITY, tid, ctx->reactivity) < 0
        || h5reader_field(ctx->reader, OUT_COVERAGE, tid, ctx->coverage) < 0) {
        snprintf(error, error_len, "%s: %s", ctx->cfg->input_path,
                 h5reader_error(ctx->reader));
        return -1;
    }

    return 0;
}

/* Writes one reference's row, with the header before the first of them, so that a run
 * scoring nothing writes nothing. */
static void print_result(const char *name, const result *r, size_t written)
{
    if (written == 0) {
        printf("reference\tpaired\tunpaired\tauroc\tauprc\tmean_paired\tmean_unpaired\n");
    }

    printf("%s\t%zu\t%zu\t%.5f\t%.5f\t%.6g\t%.6g\n", name, r->paired, r->unpaired,
           r->auroc, r->auprc, r->mean_paired, r->mean_unpaired);
}

/* ------------------------------------------------------------------------ */
/* The run                                                                   */
/* ------------------------------------------------------------------------ */

static int allocate(context *ctx, size_t cap, char *error, size_t error_len)
{
    ctx->cap        = cap;
    ctx->reactivity = malloc(out_values(OUT_REACTIVITY, cap, cap) * sizeof *ctx->reactivity);
    ctx->coverage   = malloc(out_values(OUT_COVERAGE, cap, cap) * sizeof *ctx->coverage);
    ctx->points     = malloc(cap * sizeof *ctx->points);

    if (!ctx->reactivity || !ctx->coverage || !ctx->points) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    return 0;
}

static void release(context *ctx)
{
    free(ctx->reactivity);
    free(ctx->coverage);
    free(ctx->points);
}

/* Walks the FASTA, scoring each reference that a structure is held for. */
static int score_all(context *ctx, const structures *set, char *error, size_t error_len)
{
    const char      *why    = NULL;
    cm_fasta_reader *reader = cm_fasta_open(ctx->cfg->fasta_path, &why);
    cm_fasta_record  ref;
    size_t           scored = 0;
    int32_t          tid    = 0;

    if (!reader) {
        snprintf(error, error_len, "%s: %s", ctx->cfg->fasta_path, why);
        return -1;
    }

    while (cm_fasta_next(reader, &ref) == CM_ITER_OK) {
        const structure *known = structure_named(set, ref.name);
        size_t           unpaired;
        size_t           n;
        result           one;

        if (tid >= h5reader_refs(ctx->reader)) {
            snprintf(error, error_len, "%s holds %d references, fewer than %s",
                     ctx->cfg->input_path, h5reader_refs(ctx->reader),
                     ctx->cfg->fasta_path);
            cm_fasta_close(reader);
            return -1;
        }

        /* the rows are as wide as the longest reference the file was written for, so a
         * longer one means this is not the FASTA it was counted against */
        if (ref.len > ctx->cap) {
            snprintf(error, error_len, "%s is %zu long, wider than the rows of %s",
                     ref.name, ref.len, ctx->cfg->input_path);
            cm_fasta_close(reader);
            return -1;
        }

        if (!known) {
            tid++;
            continue;
        }

        if (known->seq) {
            long at = disagreement(known->seq, &ref);

            if (at >= 0) {
                snprintf(error, error_len,
                         "%s: the structure holds %c at position %ld, where the "
                         "reference holds %c",
                         ref.name, known->seq[at], at + 1, ref.seq[at]);
                cm_fasta_close(reader);
                return -1;
            }
        }

        if (read_row(ctx, tid, error, error_len) < 0) {
            cm_fasta_close(reader);
            return -1;
        }

        n = collect(ctx, &ref, known, &unpaired);
        tid++;

        /* a reference holding one class alone gives no ranking to measure */
        if (unpaired == 0 || unpaired == n) {
            continue;
        }

        one = measure(ctx->points, n, unpaired);
        print_result(ref.name, &one, scored);
        scored++;
    }

    if (cm_fasta_error(reader)) {
        snprintf(error, error_len, "%s: %s", ctx->cfg->fasta_path,
                 cm_fasta_error(reader));
        cm_fasta_close(reader);
        return -1;
    }

    cm_fasta_close(reader);

    if (scored == 0) {
        snprintf(error, error_len, "no reference could be scored");
        return -1;
    }

    return 0;
}

int score_run(const score_config *cfg, char *error, size_t error_len)
{
    context    ctx = { .cfg = cfg };
    structures set = { 0 };
    int        status;

    if (structures_read(&set, cfg->structures_path, error, error_len) < 0) {
        structures_free(&set);
        return -1;
    }

    ctx.reader = h5reader_open(cfg->input_path);

    if (!ctx.reader || h5reader_error(ctx.reader)) {
        structures_free(&set);
        return h5reader_fail(ctx.reader, cfg->input_path, error, error_len);
    }

    if (allocate(&ctx, h5reader_capacity(ctx.reader), error, error_len) < 0) {
        status = -1;
    } else {
        status = score_all(&ctx, &set, error, error_len);
    }

    release(&ctx);
    h5reader_close(ctx.reader);
    structures_free(&set);

    return status;
}
