/* normalize.c -- one norm, divided out of every input.
 *
 * A scheme that computes the norm from the rates runs in two passes. The first reads the
 * rates every input holds and pools the values the scheme draws on. The second re-reads
 * each input and writes it out divided by the norm. The value scheme is given its norm on
 * the command line, so it makes the second pass alone.
 *
 * The pool holds the aggregate rate of the pooled channels an input has: one less the
 * product of their no-event rates. The one norm then divides every channel alike.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "normalize.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "format.h"
#include "h5reader.h"
#include "h5writer.h"
#include "progress.h"

/* How far up the pool the ubr norm is taken, as a fraction. */
#define UBR_PERCENTILE 0.90

/* The band the outlier norm averages, as fractions of the pool counted from the highest
 * value down. */
#define OUTLIER_HIGHEST 0.02
#define OUTLIER_LOWEST  0.10

/* How many rates the pool holds before it first grows. */
#define POOL_INITIAL_CAPACITY 1024

/* The channels the norm is pooled from. The termination rate is excluded, because a
 * termination is not a chemical modification. The norm divides every channel, which keeps
 * the rates comparable. */
static const fmt_channel POOLED[] = {
    FMT_CHANNEL_MISMATCHES,
    FMT_CHANNEL_INSERTIONS,
    FMT_CHANNEL_DELETIONS,
};

#define N_POOLED (sizeof POOLED / sizeof *POOLED)

/* The rates the scheme uses, gathered from every input. */
typedef struct {
    float *value;
    size_t count;
    size_t capacity;
} rate_pool;

/* ------------------------------------------------------------------------ */
/* The pool                                                                  */
/* ------------------------------------------------------------------------ */

static int rate_pool_push(rate_pool *p, float value)
{
    if (p->count == p->capacity) {
        size_t wanted = p->capacity ? p->capacity * 2 : POOL_INITIAL_CAPACITY;
        float *grown  = realloc(p->value, wanted * sizeof *grown);

        if (!grown) {
            return -1;
        }

        p->value    = grown;
        p->capacity = wanted;
    }

    p->value[p->count++] = value;
    return 0;
}

static void rate_pool_free(rate_pool *p)
{
    free(p->value);
}

/* ------------------------------------------------------------------------ */
/* Order statistics                                                          */
/* ------------------------------------------------------------------------ */

static void swap_f32(float *a, float *b)
{
    float held = *a;

    *a = *b;
    *b = held;
}

/* Puts the median of the first, middle and last value at high, so that a pool already in
 * order does not partition into one long side and one empty one. */
static void place_pivot(float *value, ptrdiff_t low, ptrdiff_t high)
{
    ptrdiff_t mid = low + ((high - low) / 2);

    if (value[mid] < value[low]) {
        swap_f32(&value[mid], &value[low]);
    }
    if (value[high] < value[low]) {
        swap_f32(&value[high], &value[low]);
    }
    if (value[high] < value[mid]) {
        swap_f32(&value[high], &value[mid]);
    }

    swap_f32(&value[mid], &value[high]);
}

/* Puts the k-th smallest value at index k, every value before it no greater and every
 * value after it no less.
 *
 * Values equal to the pivot are gathered in the same pass and then dropped from the
 * search, which keeps the cost linear for a pool of many repeated rates. */
static void select_nth(float *value, size_t n, size_t k)
{
    ptrdiff_t low  = 0;
    ptrdiff_t high = (ptrdiff_t)n - 1;
    ptrdiff_t want = (ptrdiff_t)k;

    while (low < high) {
        ptrdiff_t below = low;
        ptrdiff_t above = high;
        ptrdiff_t at    = low;
        float     pivot;

        place_pivot(value, low, high);
        pivot = value[high];

        while (at <= above) {
            if (value[at] < pivot) {
                swap_f32(&value[at], &value[below]);
                below++;
                at++;
            } else if (value[at] > pivot) {
                swap_f32(&value[at], &value[above]);
                above--;
            } else {
                at++;
            }
        }

        if (want < below) {
            high = below - 1;
        } else if (want > above) {
            low = above + 1;
        } else {
            return;
        }
    }
}

static float smallest(const float *value, size_t n)
{
    float least = value[0];

    for (size_t i = 1; i < n; i++) {
        if (value[i] < least) {
            least = value[i];
        }
    }

    return least;
}

/* Returns the value a fraction of the way up the pool, interpolating between the two
 * order statistics on either side of it. Reorders the pool. */
static double percentile(float *value, size_t n, double fraction)
{
    double position = fraction * (double)(n - 1);
    size_t below    = (size_t)position;
    double part     = position - (double)below;
    double lower;

    select_nth(value, n, below);
    lower = (double)value[below];

    if (below + 1 >= n) {
        return lower;
    }

    /* Everything past below is no less than it, so the next order statistic up is the
     * least of what remains. */
    return lower + (part * ((double)smallest(value + below + 1, n - below - 1) - lower));
}

/* ------------------------------------------------------------------------ */
/* The schemes                                                               */
/* ------------------------------------------------------------------------ */

/* The ubr and outlier schemes reorder the pool, each reading the order statistics it
 * needs out of it. */

/* Returns the rank a fraction of the way down from the highest value, counting from zero
 * and never reaching the highest value itself. */
static size_t rank_from_top(size_t n, double fraction)
{
    long at = lround(fraction * (double)n) - 1;

    return at < 1 ? 1 : (size_t)at;
}

static double ubr_norm(rate_pool *p)
{
    return p->count ? percentile(p->value, p->count, UBR_PERCENTILE) : 1.0;
}

/* Averages the band between the two fractions. This drops the highest rates as outliers
 * and computes the norm from the rates just below them. */
static double outlier_norm(rate_pool *p)
{
    size_t lowest, highest, first, last;
    double total = 0.0;

    if (p->count < 2) {
        return 1.0;
    }

    lowest  = rank_from_top(p->count, OUTLIER_LOWEST);
    highest = rank_from_top(p->count, OUTLIER_HIGHEST);

    /* Ranks counted from the top, as indices into the pool ordered upwards. */
    first = p->count - 1 - lowest;
    last  = p->count - 1 - highest;

    select_nth(p->value, p->count, first);
    select_nth(p->value + first, p->count - first, last - first);

    for (size_t i = first; i <= last; i++) {
        total += (double)p->value[i];
    }

    return total / (double)((last - first) + 1);
}

/* Returns the norm the ubr or outlier scheme computes from the pool. */
static double pooled_norm(const normalize_config *cfg, rate_pool *p)
{
    return cfg->scheme == NORM_UBR ? ubr_norm(p) : outlier_norm(p);
}

/* Returns the norm, or NaN where the number given cannot be one. A norm is a divisor, so
 * it must be above zero. */
static double usable_norm(double norm)
{
    return (isnan(norm) || norm <= 0.0) ? (double)NAN : norm;
}

/* Whether the run computes the norm. The value scheme is given its norm, so it reads no
 * rates. */
static bool computes_norm(const normalize_config *cfg)
{
    return cfg->scheme != NORM_VALUE;
}

/* Returns the number the rates are divided by. Where there is no norm, the rates are
 * divided by one and do not change. */
static double divisor(double norm)
{
    return isnan(norm) ? 1.0 : norm;
}

/* ------------------------------------------------------------------------ */
/* Arithmetic                                                                */
/* ------------------------------------------------------------------------ */

static void normalize_f32(float *row, size_t n, double norm)
{
    for (size_t i = 0; i < n; i++) {
        row[i] = (float)((double)row[i] / norm);
    }
}

/* Whether the norm divides this field. The norm divides every rate and every error.
 * Every count is left as it is. */
static bool is_normalized(fmt_field_id id)
{
    return fmt_is_channel(id);
}

/* ------------------------------------------------------------------------ */
/* Failures                                                                  */
/* ------------------------------------------------------------------------ */

static int fail_memory(char *error, size_t error_len)
{
    snprintf(error, error_len, "out of memory");
    return -1;
}

/* ------------------------------------------------------------------------ */
/* Gathering the pool                                                        */
/* ------------------------------------------------------------------------ */

/* Which of the pooled rates one input holds, and a row per rate to read into. */
typedef struct {
    fmt_field_id field[N_POOLED];
    float       *row[N_POOLED];
    size_t       n;
} pooled_rates;

/* Names the pooled rates this input holds. */
static void pooled_rates_of(pooled_rates *held, const h5reader *in)
{
    held->n = 0;

    for (size_t c = 0; c < N_POOLED; c++) {
        fmt_field_id id = FMT_CHANNEL_RATES[POOLED[c]];

        if (h5reader_holds(in, id)) {
            held->field[held->n++] = id;
        }
    }
}

/* Allocates a row for each named rate. Returns -1 where a row cannot be allocated. */
static int pooled_rates_alloc(pooled_rates *held, size_t values)
{
    for (size_t c = 0; c < held->n; c++) {
        held->row[c] = calloc(values, sizeof *held->row[c]);

        if (!held->row[c]) {
            return -1;
        }
    }

    return 0;
}

static void pooled_rates_free(pooled_rates *held)
{
    for (size_t c = 0; c < held->n; c++) {
        free(held->row[c]);
    }
}

/* Returns the rate of an event of any kind: one less the product of the no-event rates.
 * NaN in any rate makes the result NaN, so a position has an aggregate only where every
 * rate is present. */
static float aggregate_rate(const pooled_rates *held, size_t i)
{
    float none = 1.0F;

    for (size_t c = 0; c < held->n; c++) {
        none *= 1.0F - held->row[c][i];
    }

    return 1.0F - none;
}

/* Whether a position's rate joins the pool. A rate computed from few reads is
 * unreliable. A position whose coverage does not clear the floor sets no norm. */
static bool joins_pool(const normalize_config *cfg, float rate, float coverage)
{
    return isfinite(rate) && (double)coverage > cfg->min_coverage;
}

static int gather_reference(const normalize_config *cfg, rate_pool *p,
                            const pooled_rates *held, const float *coverage, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        float combined = aggregate_rate(held, i);

        if (joins_pool(cfg, combined, coverage[i]) && rate_pool_push(p, combined) < 0) {
            return -1;
        }
    }

    return 0;
}

/* Reads one reference's coverage and held rates into the rows they occupy. */
static int read_channels(h5reader *in, int32_t tid, const pooled_rates *held,
                         float *cover)
{
    if (h5reader_field(in, FMT_COVERAGE, tid, cover) < 0) {
        return -1;
    }

    for (size_t c = 0; c < held->n; c++) {
        if (h5reader_field(in, held->field[c], tid, held->row[c]) < 0) {
            return -1;
        }
    }

    return 0;
}

static int gather_input(const normalize_config *cfg, rate_pool *p, h5reader *in,
                        const char *path, char *error, size_t error_len)
{
    size_t       cap    = h5reader_capacity(in);
    size_t       values = fmt_values(FMT_MISMATCH_RATE, cap, cap);
    pooled_rates held   = { 0 };
    float       *cover  = calloc(values, sizeof *cover);
    int          status = -1;

    pooled_rates_of(&held, in);

    if (!cover || pooled_rates_alloc(&held, values) < 0) {
        fail_memory(error, error_len);
        goto done;
    }

    for (int32_t tid = 0; tid < h5reader_refs(in); tid++) {
        if (read_channels(in, tid, &held, cover) < 0) {
            h5reader_fail(in, path, error, error_len);
            goto done;
        }

        if (gather_reference(cfg, p, &held, cover, values) < 0) {
            fail_memory(error, error_len);
            goto done;
        }
    }

    status = 0;

done:
    pooled_rates_free(&held);
    free(cover);
    return status;
}

/* Reads every input in turn, so a run with a file it cannot read stops before any output
 * is created. */
static int gather(const normalize_config *cfg, const fmt_manifest *writes, rate_pool *p,
                  progress *bar, char *error, size_t error_len)
{
    fmt_request requests[FMT_N_FIELDS];
    size_t      n = fmt_requests_of(writes, requests);

    for (size_t i = 0; i < cfg->n_files; i++) {
        h5reader *in     = h5reader_open(cfg->inputs[i], requests, n);
        int       status = -1;

        if (!in) {
            return fail_memory(error, error_len);
        }

        if (h5reader_error(in)) {
            h5reader_fail(in, cfg->inputs[i], error, error_len);
        } else {
            const char *ignored = h5reader_ignored(in);

            if (ignored) {
                fprintf(stderr, "%s: ignoring %s\n", cfg->inputs[i], ignored);
            }

            status = gather_input(cfg, p, in, cfg->inputs[i], error, error_len);
        }

        h5reader_close(in);

        if (status < 0) {
            return -1;
        }

        progress_follow(bar, i + 1);
    }

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Writing one output                                                        */
/* ------------------------------------------------------------------------ */

typedef struct {
    const fmt_manifest     *manifest;               /* what it writes */
    fmt_request             requests[FMT_N_FIELDS]; /* derived from the manifest */
    size_t                  n_requests;
    double                  norm;
    bool                    writes[FMT_N_FIELDS];   /* the fields this run writes */

    h5reader   *in;
    h5writer   *out;
    const char *in_path;
    const char *out_path;
    const char *program;

    void   *row;
    int32_t n_refs;
    size_t  ref_cap;
} transfer;

static void transfer_row(const transfer *t, fmt_field_id id, size_t n)
{
    if (is_normalized(id)) {
        normalize_f32(t->row, n, divisor(t->norm));
    }
}

static int transfer_reference(const transfer *t, int32_t tid, char *error,
                              size_t error_len)
{
    for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
        if (!FMT_FIELDS[id].per_ref || !fmt_wanted(id, t->writes)) {
            continue;
        }

        if (h5reader_field(t->in, id, tid, t->row) < 0) {
            return h5reader_fail(t->in, t->in_path, error, error_len);
        }

        transfer_row(t, id, fmt_values(id, t->ref_cap, t->ref_cap));

        if (h5writer_row(t->out, id, tid, t->row) < 0) {
            return h5writer_fail(t->out, t->out_path, error, error_len);
        }
    }

    return 0;
}

static int transfer_totals(const transfer *t, char *error, size_t error_len)
{
    for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
        size_t total;

        if (FMT_FIELDS[id].per_ref || !h5reader_holds(t->in, id)) {
            continue;
        }

        if (h5reader_total(t->in, id, &total) < 0) {
            return h5reader_fail(t->in, t->in_path, error, error_len);
        }

        if (h5writer_total(t->out, id, total) < 0) {
            return h5writer_fail(t->out, t->out_path, error, error_len);
        }
    }

    return 0;
}

static int transfer_file(const transfer *t, char *error, size_t error_len)
{
    for (int32_t tid = 0; tid < t->n_refs; tid++) {
        if (transfer_reference(t, tid, error, error_len) < 0) {
            return -1;
        }
    }

    if (transfer_totals(t, error, error_len) < 0) {
        return -1;
    }

    if (h5writer_value(t->out, FMT_NORM, t->norm) < 0) {
        return h5writer_fail(t->out, t->out_path, error, error_len);
    }

    return 0;
}

/* Clears from writes every field that depends on one the input does not have. What
 * remains is read, copied and written alike. */
static void drop_absent_fields(transfer *t)
{
    for (size_t i = 0; i < t->manifest->n_fields; i++) {
        const fmt_written *field = &t->manifest->fields[i];

        for (const fmt_field_id *dep = field->depends;
             dep && *dep != FMT_N_FIELDS; dep++) {
            if (!h5reader_holds(t->in, *dep)) {
                t->writes[field->id] = false;
            }
        }
    }
}

static int open_transfer(transfer *t, bool may_replace, char *error, size_t error_len)
{
    t->n_requests = fmt_requests_of(t->manifest, t->requests);

    t->in = h5reader_open(t->in_path, t->requests, t->n_requests);

    if (!t->in) {
        return fail_memory(error, error_len);
    }

    if (h5reader_error(t->in)) {
        return h5reader_fail(t->in, t->in_path, error, error_len);
    }

    drop_absent_fields(t);

    t->n_refs  = h5reader_refs(t->in);
    t->ref_cap = h5reader_capacity(t->in);
    t->row     = calloc(fmt_widest(t->ref_cap, t->writes), fmt_widest_bytes());

    if (!t->row) {
        return fail_memory(error, error_len);
    }

    t->out = h5writer_create(t->out_path, t->program, t->n_refs, t->ref_cap,
                             may_replace, t->writes, false);
    if (!t->out) {
        return fail_memory(error, error_len);
    }

    return h5writer_error(t->out) ? h5writer_fail(t->out, t->out_path, error, error_len) : 0;
}

static void transfer_teardown(transfer *t)
{
    h5writer_close(t->out);
    h5reader_close(t->in);
    free(t->row);
}

static int write_output(const normalize_config *cfg, size_t which, const char *program,
                        const fmt_manifest *writes,
                        double norm, char *error, size_t error_len)
{
    transfer t = {
        .manifest = writes,
        .norm   = norm,
        .in_path  = cfg->inputs[which],
        .out_path = cfg->outputs[which],
        .program  = program,
    };
    bool may_replace = false;
    int  status      = -1;

    fmt_selection(writes, t.writes);

    /* Checked again here rather than reused from check_outputs, so that a file created
     * at the path since then is found. */
    if (h5writer_may_replace(t.out_path, cfg->overwrite, &may_replace, error,
                             error_len) < 0) {
        return -1;
    }

    if (open_transfer(&t, may_replace, error, error_len) == 0) {
        status = transfer_file(&t, error, error_len);
    }

    transfer_teardown(&t);
    return status;
}

/* ------------------------------------------------------------------------ */
/* Assembly                                                                  */
/* ------------------------------------------------------------------------ */

/* Refuses every output path the run could not write, before the first is created, so that
 * a bad path late in the list does not waste the earlier outputs. Whether each may be
 * replaced is checked again when the file is created. */
static int check_outputs(const normalize_config *cfg, char *error, size_t error_len)
{
    for (size_t i = 0; i < cfg->n_files; i++) {
        bool may_replace = false;

        if (h5writer_may_replace(cfg->outputs[i], cfg->overwrite, &may_replace, error,
                                 error_len) < 0) {
            return -1;
        }
    }

    return 0;
}

/* How many passes over the inputs a run makes. */
static uint64_t passes(const normalize_config *cfg)
{
    return computes_norm(cfg) ? 2 : 1;
}

/* How many units of progress the run counts before it writes the outputs. */
static uint64_t units_before_writing(const normalize_config *cfg)
{
    return (passes(cfg) - 1) * (uint64_t)cfg->n_files;
}

static int write_outputs(const normalize_config *cfg, const char *program,
                         const fmt_manifest *writes, double norm,
                         progress *bar, char *error, size_t error_len)
{
    uint64_t done = units_before_writing(cfg);

    for (size_t i = 0; i < cfg->n_files; i++) {
        if (write_output(cfg, i, program, writes, norm, error, error_len) < 0) {
            return -1;
        }

        progress_follow(bar, done + i + 1);
    }

    return 0;
}

/* Writes the norm every input is divided by to norm. The value scheme reads that norm
 * from the command line. Every other scheme computes it from the rates, which is the
 * first pass. Returns 0, or -1 with a description in error. */
static int norm_of(const normalize_config *cfg, const fmt_manifest *writes, progress *bar,
                   double *norm, char *error, size_t error_len)
{
    rate_pool p = { 0 };
    int       status;

    if (!computes_norm(cfg)) {
        *norm = usable_norm(cfg->value);
        return 0;
    }

    status = gather(cfg, writes, &p, bar, error, error_len);

    if (status == 0) {
        *norm = usable_norm(pooled_norm(cfg, &p));
    }

    rate_pool_free(&p);
    return status;
}

int normalize_run(const normalize_config *cfg, const char *program,
                  const fmt_manifest *writes, char *error,
                  size_t error_len)
{
    int       status = -1;
    double    norm   = 0.0;
    progress *bar;

    if (check_outputs(cfg, error, error_len) < 0) {
        return -1;
    }

    /* One unit per input per pass. */
    bar = progress_start(passes(cfg) * (uint64_t)cfg->n_files);

    if (norm_of(cfg, writes, bar, &norm, error, error_len) == 0) {
        status = write_outputs(cfg, program, writes, norm, bar, error, error_len);
    }

    progress_finish(bar);
    return status;
}
