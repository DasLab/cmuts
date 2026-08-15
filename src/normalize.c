/* normalize.c -- one scale, pooled over every input, divided out of each.
 *
 * Done in two passes. The first reads the rates every input holds and pools the values
 * the scheme draws on; the second re-reads each input and writes it out divided by the
 * scale.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "normalize.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "h5reader.h"
#include "h5writer.h"
#include "output.h"

/* The name written into the output as the program that produced it. */
#define NORMALIZE_PROGRAM "cmuts-norm"

/* The rate the ubr scale sits at, as a fraction of the way up the pool. */
#define UBR_PERCENTILE 0.90

/* The band the outlier scale averages, as fractions of the pool counted from the highest
 * value down. */
#define OUTLIER_HIGHEST 0.02
#define OUTLIER_LOWEST  0.10

/* Rates the pool holds before it first grows. */
#define POOL_INITIAL_CAPACITY 1024

/* The rates the scheme draws on, gathered from every input. */
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
 * search, which is what holds a pool of many repeated rates to a linear cost. */
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

/* Returns the value a fraction of the way up the pool, interpolating between the two order
 * statistics it falls between. Reorders the pool. */
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

/* Both schemes reorder the pool, each reading the order statistics it needs out of it. */

/* Returns the rank a fraction of the way down from the highest value, counting from zero
 * and never reaching the highest value itself. */
static size_t rank_from_top(size_t n, double fraction)
{
    long at = lround(fraction * (double)n) - 1;

    return at < 1 ? 1 : (size_t)at;
}

static double ubr_factor(rate_pool *p)
{
    return p->count ? percentile(p->value, p->count, UBR_PERCENTILE) : 1.0;
}

/* Averages the band between the two fractions, which drops the highest rates as outliers
 * and takes the scale from what sits just below them. */
static double outlier_factor(rate_pool *p)
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

static double pooled_factor(const normalize_config *cfg, rate_pool *p)
{
    double factor = cfg->scheme == NORM_UBR ? ubr_factor(p) : outlier_factor(p);

    /* A factor that is not above zero leaves the rates unscaled. */
    return (isnan(factor) || factor <= 0.0) ? 1.0 : factor;
}

/* ------------------------------------------------------------------------ */
/* Arithmetic                                                                */
/* ------------------------------------------------------------------------ */

static void scale_f32(float *row, size_t n, double factor)
{
    for (size_t i = 0; i < n; i++) {
        row[i] = (float)((double)row[i] / factor);
    }
}

/* Holds every value within the bounds, leaving NaN as it is and leaving a bound that is
 * NaN unapplied. */
static void clip_f32(float *row, size_t n, double below, double above)
{
    bool  has_below = !isnan(below);
    bool  has_above = !isnan(above);
    float low       = (float)below;
    float high      = (float)above;

    for (size_t i = 0; i < n; i++) {
        if (has_below && row[i] < low) {
            row[i] = low;
        }
        if (has_above && row[i] > high) {
            row[i] = high;
        }
    }
}

/* Whether the scale divides this field. The rate and its error take it; every count is
 * left as it stands. */
static bool is_scaled(out_field_id id)
{
    return id == OUT_REACTIVITY || id == OUT_ERROR;
}

/* ------------------------------------------------------------------------ */
/* Failures                                                                  */
/* ------------------------------------------------------------------------ */

static int fail_input(const h5reader *in, const char *path, char *error,
                      size_t error_len)
{
    const char *why = h5reader_error(in);

    snprintf(error, error_len, "%s: %s", path, why ? why : "unable to read it");
    return -1;
}

static int fail_output(const h5writer *out, const char *path, char *error,
                       size_t error_len)
{
    const char *why = h5writer_error(out);

    snprintf(error, error_len, "%s: %s", path, why ? why : "unable to write it");
    return -1;
}

static int fail_memory(char *error, size_t error_len)
{
    snprintf(error, error_len, "out of memory");
    return -1;
}

/* ------------------------------------------------------------------------ */
/* Gathering the pool                                                        */
/* ------------------------------------------------------------------------ */

/* Whether a position's rate joins the pool. ubr takes only the positions whose coverage
 * clears the floor; outlier takes every rate there is. */
static bool joins_pool(const normalize_config *cfg, float rate, float coverage)
{
    if (!isfinite(rate)) {
        return false;
    }

    return cfg->scheme != NORM_UBR || (double)coverage > cfg->min_coverage;
}

static int gather_reference(const normalize_config *cfg, rate_pool *p, const float *rate,
                            const float *coverage, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        if (joins_pool(cfg, rate[i], coverage[i]) && rate_pool_push(p, rate[i]) < 0) {
            return -1;
        }
    }

    return 0;
}

static int gather_input(const normalize_config *cfg, rate_pool *p, h5reader *in,
                        const char *path, char *error, size_t error_len)
{
    size_t  cap    = h5reader_capacity(in);
    size_t  values = out_values(OUT_REACTIVITY, cap, cap);
    float  *rate   = calloc(values, sizeof *rate);
    float  *cover  = calloc(values, sizeof *cover);
    int     status = -1;

    if (!rate || !cover) {
        fail_memory(error, error_len);
        goto done;
    }

    for (int32_t tid = 0; tid < h5reader_refs(in); tid++) {
        if (h5reader_field(in, OUT_REACTIVITY, tid, rate) < 0 ||
            h5reader_field(in, OUT_COVERAGE, tid, cover) < 0) {
            fail_input(in, path, error, error_len);
            goto done;
        }

        if (gather_reference(cfg, p, rate, cover, values) < 0) {
            fail_memory(error, error_len);
            goto done;
        }
    }

    status = 0;

done:
    free(rate);
    free(cover);
    return status;
}

/* Reads every input in turn, so a file that cannot be read fails before any output is
 * created. */
static int gather(const normalize_config *cfg, rate_pool *p, char *error,
                  size_t error_len)
{
    for (size_t i = 0; i < cfg->n_files; i++) {
        h5reader *in     = h5reader_open(cfg->inputs[i]);
        int       status = -1;

        if (!in) {
            return fail_memory(error, error_len);
        }

        if (h5reader_error(in)) {
            fail_input(in, cfg->inputs[i], error, error_len);
        } else {
            status = gather_input(cfg, p, in, cfg->inputs[i], error, error_len);
        }

        h5reader_close(in);

        if (status < 0) {
            return -1;
        }
    }

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Writing one output                                                        */
/* ------------------------------------------------------------------------ */

typedef struct {
    const normalize_config *cfg;
    double                  factor;

    h5reader   *in;
    h5writer   *out;
    const char *in_path;
    const char *out_path;

    void   *row;
    int32_t n_refs;
    size_t  ref_cap;
} transfer;

static void transfer_row(const transfer *t, out_field_id id, size_t n)
{
    if (!is_scaled(id)) {
        return;
    }

    scale_f32(t->row, n, t->factor);

    if (id == OUT_REACTIVITY) {
        clip_f32(t->row, n, t->cfg->clip_below, t->cfg->clip_above);
    }
}

static int transfer_reference(const transfer *t, int32_t tid, char *error,
                              size_t error_len)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (!OUT_FIELDS[id].per_ref) {
            continue;
        }

        if (h5reader_field(t->in, id, tid, t->row) < 0) {
            return fail_input(t->in, t->in_path, error, error_len);
        }

        transfer_row(t, id, out_values(id, t->ref_cap, t->ref_cap));

        if (h5writer_row(t->out, id, tid, t->row) < 0) {
            return fail_output(t->out, t->out_path, error, error_len);
        }
    }

    return 0;
}

static int transfer_totals(const transfer *t, char *error, size_t error_len)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        size_t total;

        if (OUT_FIELDS[id].per_ref) {
            continue;
        }

        if (h5reader_total(t->in, id, &total) < 0) {
            return fail_input(t->in, t->in_path, error, error_len);
        }

        if (h5writer_total(t->out, id, total) < 0) {
            return fail_output(t->out, t->out_path, error, error_len);
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

    if (h5writer_attribute(t->out, NORMALIZE_ATTRIBUTE, t->factor) < 0) {
        return fail_output(t->out, t->out_path, error, error_len);
    }

    return 0;
}

static int open_transfer(transfer *t, bool may_replace, char *error, size_t error_len)
{
    t->in = h5reader_open(t->in_path);

    if (!t->in) {
        return fail_memory(error, error_len);
    }

    if (h5reader_error(t->in)) {
        return fail_input(t->in, t->in_path, error, error_len);
    }

    t->n_refs  = h5reader_refs(t->in);
    t->ref_cap = h5reader_capacity(t->in);
    t->row     = calloc(out_widest(t->ref_cap), out_widest_bytes());

    if (!t->row) {
        return fail_memory(error, error_len);
    }

    t->out = h5writer_create(t->out_path, NORMALIZE_PROGRAM, t->n_refs, t->ref_cap,
                             may_replace);
    if (!t->out) {
        return fail_memory(error, error_len);
    }

    return h5writer_error(t->out) ? fail_output(t->out, t->out_path, error, error_len) : 0;
}

static void transfer_teardown(transfer *t)
{
    h5writer_close(t->out);
    h5reader_close(t->in);
    free(t->row);
}

static int write_output(const normalize_config *cfg, size_t which, double factor,
                        char *error, size_t error_len)
{
    transfer t = {
        .cfg      = cfg,
        .factor   = factor,
        .in_path  = cfg->inputs[which],
        .out_path = cfg->outputs[which],
    };
    bool may_replace = false;
    int  status      = -1;

    /* Asked again here rather than carried over from check_outputs, so that a file
     * appearing at the path since then is seen. */
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
 * a bad path late in the list costs none of the earlier ones. Whether each may be replaced
 * is settled again at the create itself. */
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

static int write_outputs(const normalize_config *cfg, double factor, char *error,
                         size_t error_len)
{
    for (size_t i = 0; i < cfg->n_files; i++) {
        if (write_output(cfg, i, factor, error, error_len) < 0) {
            return -1;
        }
    }

    return 0;
}

int normalize_run(const normalize_config *cfg, char *error, size_t error_len)
{
    rate_pool p      = { 0 };
    int       status = -1;

    if (check_outputs(cfg, error, error_len) < 0) {
        return -1;
    }

    if (gather(cfg, &p, error, error_len) == 0) {
        status = write_outputs(cfg, pooled_factor(cfg, &p), error, error_len);
    }

    rate_pool_free(&p);
    return status;
}
