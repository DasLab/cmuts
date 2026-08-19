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

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
#include "fasta.h"
#include "h5reader.h"
#include "nuc.h"
#include "output.h"

/* What a dot bracket record may hold. A base is paired under any bracket, open under a
 * dot, and unknown under anything else, which is left out of the scoring. */
#define BRACKETS "()[]{}<>"
#define UNPAIRED '.'

/* The longest line a record may hold, which bounds the longest structure that can be
 * read, and the longest name kept from a header. */
#define LINE_MAX_LEN 4096
#define NAME_MAX_LEN 256

/* Structures the set holds before it first grows. */
#define STRUCTURES_INITIAL_CAPACITY 64

/* ------------------------------------------------------------------------ */
/* The structures                                                            */
/* ------------------------------------------------------------------------ */

/* One reference's structure. The name and the pairing own their memory. */
typedef struct {
    char  *name;
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
    size_t     cap = set->cap ? set->cap * 2 : STRUCTURES_INITIAL_CAPACITY;
    structure *at  = realloc(set->at, cap * sizeof *at);

    if (!at) {
        return -1;
    }

    set->at  = at;
    set->cap = cap;

    return 0;
}

/* Returns whether the line is a pairing and not a sequence. Every character of a
 * pairing is a bracket, a dot or a gap, and a sequence carries none of them. */
static bool is_pairing(const char *line)
{
    size_t len = strlen(line);

    return len > 0 && strspn(line, BRACKETS ".-") == len;
}

/* Cuts the line at its first line ending. */
static void trim(char *line)
{
    char *end = strpbrk(line, "\r\n");

    if (end) {
        *end = '\0';
    }
}

/* Keeps one record, copying it, since the lines it came from are reused. */
static int structures_add(structures *set, const char *name, const char *pairing)
{
    structure *at;

    if (set->n == set->cap && structures_grow(set) < 0) {
        return -1;
    }

    at          = &set->at[set->n];
    at->name    = strdup(name);
    at->pairing = strdup(pairing);
    at->len     = strlen(pairing);

    if (!at->name || !at->pairing) {
        free(at->name);
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
    char  line[LINE_MAX_LEN];
    char  name[NAME_MAX_LEN] = { 0 };
    FILE *file               = fopen(path, "r");
    int   number             = 0;

    if (!file) {
        snprintf(error, error_len, "%s: cannot be opened", path);
        return -1;
    }

    while (fgets(line, sizeof line, file)) {
        number++;

        if (!strchr(line, '\n') && !feof(file)) {
            snprintf(error, error_len, "%s:%d: the line is too long", path, number);
            fclose(file);
            return -1;
        }

        trim(line);

        if (line[0] == '>') {
            snprintf(name, sizeof name, "%s", line + 1 + strspn(line + 1, " \t"));
            name[strcspn(name, " \t")] = '\0';
            continue;
        }

        /* a record may carry its sequence before the pairing, which is passed over */
        if (line[0] == '\0' || !is_pairing(line)) {
            continue;
        }

        if (name[0] == '\0') {
            snprintf(error, error_len, "%s:%d: a pairing before any name", path, number);
            fclose(file);
            return -1;
        }

        if (structures_add(set, name, line) < 0) {
            snprintf(error, error_len, "%s: out of memory", path);
            fclose(file);
            return -1;
        }

        name[0] = '\0';
    }

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
    double found    = 0.0;   /* unpaired bases at or above the value */
    double seen     = 0.0;   /* bases at or above it */
    double previous = 0.0;   /* recall before this run */
    double sum      = 0.0;
    size_t i        = n;

    while (i > 0) {
        size_t j = i;
        double recall;

        /* a whole run of one value is taken at once, so the order within it, which the
         * sort does not fix, cannot reach the result */
        while (j > 0 && at[j - 1].value == at[i - 1].value) {
            found += at[j - 1].open ? 1.0 : 0.0;
            seen  += 1.0;
            j--;
        }

        recall = found / (double)unpaired;
        sum   += (recall - previous) * (found / seen);
        previous = recall;
        i = j;
    }

    return sum;
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
    FILE               *out;
    h5reader           *reader;
    float              *reactivity;
    float              *coverage;
    point              *points;
    size_t              cap;
} context;

/* Returns the bit standing for a base, or 0 for a base the alphabet does not name, such
 * as an ambiguous one. T and U give the same bit. */
static int bit_of(char base)
{
    switch (nuc_from_char(base)) {
        case NUC_A: return SCORE_BASE_A;
        case NUC_C: return SCORE_BASE_C;
        case NUC_G: return SCORE_BASE_G;
        case NUC_T: return SCORE_BASE_U;
        default:    return 0;
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

    /* the structure and the reference are the same length, which score_all checks */
    for (size_t i = 0; i < ref->len; i++) {
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

/* Names the columns. Written before the first row, so that a run scoring nothing
 * writes nothing. */
static void print_header(FILE *out)
{
    fprintf(out, "reference,paired,unpaired,auroc,auprc,mean_paired,mean_unpaired\n");
}

static void print_result(FILE *out, const char *name, const result *r)
{
    fprintf(out, "%s,%zu,%zu,%.5f,%.5f,%.6g,%.6g\n", name, r->paired, r->unpaired,
            r->auroc, r->auprc, r->mean_paired, r->mean_unpaired);
}

/* ------------------------------------------------------------------------ */
/* The run                                                                   */
/* ------------------------------------------------------------------------ */

static int fail_memory(char *error, size_t error_len)
{
    snprintf(error, error_len, "out of memory");
    return -1;
}

static int allocate(context *ctx, size_t cap, char *error, size_t error_len)
{
    ctx->cap        = cap;
    ctx->reactivity = malloc(out_values(OUT_REACTIVITY, cap, cap) * sizeof *ctx->reactivity);
    ctx->coverage   = malloc(out_values(OUT_COVERAGE, cap, cap) * sizeof *ctx->coverage);
    ctx->points     = malloc(cap * sizeof *ctx->points);

    if (!ctx->reactivity || !ctx->coverage || !ctx->points) {
        return fail_memory(error, error_len);
    }

    return 0;
}

static void release(context *ctx)
{
    free(ctx->reactivity);
    free(ctx->coverage);
    free(ctx->points);
}

/* Refuses a run that scored nothing, naming the lengths where they are the cause. */
static int check_scored(size_t scored, size_t skipped, char *error, size_t error_len)
{
    if (scored > 0) {
        return 0;
    }

    if (skipped > 0) {
        snprintf(error, error_len,
                 "no reference could be scored; %zu structures are the wrong length",
                 skipped);
        return -1;
    }

    snprintf(error, error_len, "no reference could be scored");

    return -1;
}

/* Walks the FASTA, scoring each reference that a structure is held for. */
static int score_all(context *ctx, const structures *set, char *error, size_t error_len)
{
    const char      *why    = NULL;
    cm_fasta_reader *reader = cm_fasta_open(ctx->cfg->fasta_path, &why);
    cm_fasta_record  ref;
    size_t           scored  = 0;
    size_t           skipped = 0;   /* structures whose length is not the reference's */
    int32_t          tid     = 0;
    int              status  = -1;

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
            goto done;
        }

        if (!known) {
            tid++;
            continue;
        }

        /* a structure of another length describes another molecule, so scoring the part
         * that overlaps would give a number that means nothing */
        if (known->len != ref.len) {
            fprintf(stderr, "%s: the structure is %zu long and the reference %zu; "
                            "not scored\n", ref.name, known->len, ref.len);
            skipped++;
            tid++;
            continue;
        }

        /* the rows are as wide as the longest reference the file was written for, so a
         * longer one means this is not the FASTA it was counted against */
        if (ref.len > ctx->cap) {
            snprintf(error, error_len, "%s is %zu long, wider than the rows of %s",
                     ref.name, ref.len, ctx->cfg->input_path);
            goto done;
        }

        if (read_row(ctx, tid, error, error_len) < 0) {
            goto done;
        }

        n = collect(ctx, &ref, known, &unpaired);
        tid++;

        /* a reference holding one class alone gives no ranking to measure */
        if (unpaired == 0 || unpaired == n) {
            continue;
        }

        one = measure(ctx->points, n, unpaired);

        if (scored == 0) {
            print_header(ctx->out);
        }

        print_result(ctx->out, ref.name, &one);
        scored++;
    }

    if (cm_fasta_error(reader)) {
        snprintf(error, error_len, "%s: %s", ctx->cfg->fasta_path,
                 cm_fasta_error(reader));
        goto done;
    }

    status = check_scored(scored, skipped, error, error_len);

done:
    cm_fasta_close(reader);

    return status;
}

int score_run(const score_config *cfg, FILE *out, char *error, size_t error_len)
{
    context    ctx = { .cfg = cfg, .out = out };
    structures set = { 0 };
    int        status;

    if (structures_read(&set, cfg->structures_path, error, error_len) < 0) {
        structures_free(&set);
        return -1;
    }

    ctx.reader = h5reader_open(cfg->input_path);

    if (!ctx.reader) {
        structures_free(&set);
        return fail_memory(error, error_len);
    }

    if (h5reader_error(ctx.reader)) {
        status = h5reader_fail(ctx.reader, cfg->input_path, error, error_len);
        h5reader_close(ctx.reader);
        structures_free(&set);
        return status;
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
