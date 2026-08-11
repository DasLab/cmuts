/* subtract.c -- combining two outputs into one, field by field.
 *
 * Both inputs are read a row at a time and the result written the same way, so
 * memory is bounded by the longest reference and not by the size of the files.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "subtract.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "h5reader.h"
#include "h5writer.h"
#include "output.h"

typedef struct {
    const subtract_config *cfg;

    h5reader *treated;
    h5reader *untreated;
    h5writer *out;

    double *left;      /* one row of the treated file */
    double *right;     /* the same row of the untreated file */
    double *result;    /* the two combined */

    int32_t n_refs;
    size_t  ref_cap;
    size_t  totals[OUT_N_FIELDS];  /* the fields with no row, combined on opening */
} subtraction;

/* ------------------------------------------------------------------------ */
/* Failures                                                                  */
/* ------------------------------------------------------------------------ */

static int fail_input(char *error, size_t error_len, const char *path,
                      const h5reader *r)
{
    const char *why = h5reader_error(r);

    snprintf(error, error_len, "%s: %s", path, why ? why : "unable to read it");
    return -1;
}

static int fail_output(const subtraction *s, char *error, size_t error_len)
{
    const char *why = h5writer_error(s->out);

    snprintf(error, error_len, "%s: %s", s->cfg->output_path,
             why ? why : "unable to write it");
    return -1;
}

/* ------------------------------------------------------------------------ */
/* Arithmetic                                                                */
/* ------------------------------------------------------------------------ */

/* Combines one value from each input under the given rule. NaN in either input
 * carries through all three, marking the position unmeasured: a difference of rates
 * is known only where both rates are. */
static double combine(out_combine how, double treated, double untreated)
{
    switch (how) {
        case OUT_ADD:       return treated + untreated;
        case OUT_SUBTRACT:  return treated - untreated;
        case OUT_PROPAGATE: return sqrt(treated * treated + untreated * untreated);
    }

    return (double)NAN;
}

static void combine_row(out_combine how, const double *treated,
                        const double *untreated, double *out, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        out[i] = combine(how, treated[i], untreated[i]);
    }
}

/* ------------------------------------------------------------------------ */
/* Rows                                                                      */
/* ------------------------------------------------------------------------ */

static int subtract_field(subtraction *s, out_field_id id, int32_t tid,
                          char *error, size_t error_len)
{
    size_t width = out_values(id, s->ref_cap, s->ref_cap);

    if (h5reader_field(s->treated, id, tid, s->left) < 0) {
        return fail_input(error, error_len, s->cfg->treated_path, s->treated);
    }

    if (h5reader_field(s->untreated, id, tid, s->right) < 0) {
        return fail_input(error, error_len, s->cfg->untreated_path, s->untreated);
    }

    combine_row(OUT_FIELDS[id].combine, s->left, s->right, s->result, width);

    /* Written at the capacity rather than at the reference's own length, an output
     * recording no reference lengths. The columns past a reference are NaN in both
     * inputs and stay NaN through the arithmetic, which is the mark the writer would
     * otherwise apply itself. */
    if (h5writer_field(s->out, id, tid, s->ref_cap, s->result) < 0) {
        return fail_output(s, error, error_len);
    }

    return 0;
}

static int subtract_reference(subtraction *s, int32_t tid, char *error,
                              size_t error_len)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (!OUT_FIELDS[id].per_ref) {
            continue;
        }

        if (subtract_field(s, id, tid, error, error_len) < 0) {
            return -1;
        }
    }

    return 0;
}

static int subtract_references(subtraction *s, char *error, size_t error_len)
{
    for (int32_t tid = 0; tid < s->n_refs; tid++) {
        if (subtract_reference(s, tid, error, error_len) < 0) {
            return -1;
        }
    }

    return 0;
}

/* Combines the fields belonging to the run rather than to any reference, under the same
 * rules their per-reference neighbors follow.
 *
 * Read while the inputs are being opened rather than when the result is written, so that a
 * file missing one is refused before the output is created. */
static int read_totals(subtraction *s, char *error, size_t error_len)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        size_t treated, untreated;

        if (OUT_FIELDS[id].per_ref) {
            continue;
        }

        if (h5reader_total(s->treated, id, &treated) < 0) {
            return fail_input(error, error_len, s->cfg->treated_path, s->treated);
        }

        if (h5reader_total(s->untreated, id, &untreated) < 0) {
            return fail_input(error, error_len, s->cfg->untreated_path, s->untreated);
        }

        s->totals[id] = (size_t)combine(OUT_FIELDS[id].combine, (double)treated,
                                        (double)untreated);
    }

    return 0;
}

static int write_totals(subtraction *s, char *error, size_t error_len)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (OUT_FIELDS[id].per_ref) {
            continue;
        }

        if (h5writer_total(s->out, id, s->totals[id]) < 0) {
            return fail_output(s, error, error_len);
        }
    }

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Assembly                                                                  */
/* ------------------------------------------------------------------------ */

/* Checks that the two files describe the same references, and records their shape.
 * An output names no references -- a row is identified by its position -- so only
 * the shape can be compared: the same number of rows at the same capacity. */
static int check_agreement(subtraction *s, char *error, size_t error_len)
{
    if (h5reader_refs(s->treated) != h5reader_refs(s->untreated)) {
        snprintf(error, error_len,
                 "%s holds %d references and %s holds %d",
                 s->cfg->treated_path, h5reader_refs(s->treated),
                 s->cfg->untreated_path, h5reader_refs(s->untreated));
        return -1;
    }

    if (h5reader_capacity(s->treated) != h5reader_capacity(s->untreated)) {
        snprintf(error, error_len,
                 "%s is %zu bases wide and %s is %zu; the two were counted "
                 "against different references",
                 s->cfg->treated_path, h5reader_capacity(s->treated),
                 s->cfg->untreated_path, h5reader_capacity(s->untreated));
        return -1;
    }

    s->n_refs  = h5reader_refs(s->treated);
    s->ref_cap = h5reader_capacity(s->treated);

    return 0;
}

static int open_inputs(subtraction *s, char *error, size_t error_len)
{
    s->treated   = h5reader_open(s->cfg->treated_path);
    s->untreated = h5reader_open(s->cfg->untreated_path);

    if (!s->treated || !s->untreated) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    if (h5reader_error(s->treated)) {
        return fail_input(error, error_len, s->cfg->treated_path, s->treated);
    }

    if (h5reader_error(s->untreated)) {
        return fail_input(error, error_len, s->cfg->untreated_path, s->untreated);
    }

    if (check_agreement(s, error, error_len) < 0) {
        return -1;
    }

    return read_totals(s, error, error_len);
}

/* Three buffers, each one row of the widest field, so that each holds a row of any
 * field. */
static int build_buffers(subtraction *s, char *error, size_t error_len)
{
    size_t width = out_widest(s->ref_cap);

    s->left   = calloc(width, sizeof *s->left);
    s->right  = calloc(width, sizeof *s->right);
    s->result = calloc(width, sizeof *s->result);

    if (!s->left || !s->right || !s->result) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    return 0;
}

static int open_output(subtraction *s, bool may_replace, char *error,
                       size_t error_len)
{
    s->out = h5writer_create(s->cfg->output_path, s->n_refs, s->ref_cap,
                             may_replace);
    if (!s->out) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    return h5writer_error(s->out) ? fail_output(s, error, error_len) : 0;
}

static void subtraction_teardown(subtraction *s)
{
    h5writer_close(s->out);
    h5reader_close(s->untreated);
    h5reader_close(s->treated);

    free(s->result);
    free(s->right);
    free(s->left);
}

int subtract_run(const subtract_config *cfg, char *error, size_t error_len)
{
    subtraction s           = { .cfg = cfg };
    bool        may_replace = false;
    int         status      = -1;

    if (h5writer_may_replace(cfg->output_path, cfg->overwrite, &may_replace,
                             error, error_len) < 0) {
        return -1;
    }

    /* Every way an input can be wrong is found before the output is created, so a run
     * that refuses its inputs leaves whatever is at that path alone. */
    if (open_inputs(&s, error, error_len) == 0 &&
        build_buffers(&s, error, error_len) == 0 &&
        open_output(&s, may_replace, error, error_len) == 0 &&
        subtract_references(&s, error, error_len) == 0) {
        status = write_totals(&s, error, error_len);
    }

    subtraction_teardown(&s);
    return status;
}
