/* subtract.c -- combining several outputs into one, field by field.
 *
 * Every input is read a row at a time and the result written the same way, so memory is
 * bounded by the longest reference and not by the size of the files.
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
#include "version.h"

typedef enum {
    SUB_TREATED,
    SUB_UNTREATED,
    SUB_DENATURED,
    SUB_N_INPUTS,
} sub_input;

typedef enum {
    SUB_UNCONTROLLED,
    SUB_CONTROLLED,
    SUB_N_MODES,
} sub_mode;

typedef enum {
    SUB_SUM,
    SUB_DIFFERENCE,
    SUB_QUADRATURE,
    SUB_RATIO,
    SUB_RATIO_ERROR,
    SUB_N_RULES,
} sub_rule;

typedef struct {
    const char *name;
    const char *detail;  /* what the rule does to the values, in one phrase */
} sub_rule_description;

static const sub_rule_description DESCRIPTIONS[SUB_N_RULES] = {
    [SUB_SUM] = {
        .name   = "sum",
        .detail = "every input added",
    },
    [SUB_DIFFERENCE] = {
        .name   = "difference",
        .detail = "the untreated value taken from the treated one",
    },
    [SUB_QUADRATURE] = {
        .name   = "quadrature",
        .detail = "the two added in quadrature",
    },
    [SUB_RATIO] = {
        .name   = "ratio",
        .detail = "that difference over the control's value, and missing where the "
                  "control measured nothing to divide by",
    },
    [SUB_RATIO_ERROR] = {
        .name   = "ratio error",
        .detail = "the error of that ratio, taking the three runs as independent",
    },
};

static const sub_rule RULES[OUT_N_FIELDS][SUB_N_MODES] = {
    [OUT_COVERAGE]   = { SUB_SUM,        SUB_SUM         },
    [OUT_REACTIVITY] = { SUB_DIFFERENCE, SUB_RATIO       },
    [OUT_ERROR]      = { SUB_QUADRATURE, SUB_RATIO_ERROR },
    [OUT_LENGTHS]    = { SUB_SUM,        SUB_SUM         },
    [OUT_READS]      = { SUB_SUM,        SUB_SUM         },
    [OUT_REJECTED]   = { SUB_SUM,        SUB_SUM         },
    [OUT_UNMAPPED]   = { SUB_SUM,        SUB_SUM         },
};

/* One reference's rows for one file, each in the type its field is stored as. A field with
 * no row has no buffer here. */
typedef struct {
    void *field[OUT_N_FIELDS];
} row_set;

typedef struct {
    const subtract_config *cfg;

    h5reader   *input[SUB_N_INPUTS];
    const char *path[SUB_N_INPUTS];
    size_t      n_inputs;

    /* Set with n_inputs and never apart from it: a controlled rule reads the third input,
     * and nothing else keeps that read in bounds. */
    sub_mode mode;

    h5writer *out;

    row_set rows[SUB_N_INPUTS];
    void   *result;

    int32_t n_refs;
    size_t  ref_cap;
    size_t  totals[OUT_N_FIELDS];  /* only the fields with no row are filled */
} subtraction;

/* ------------------------------------------------------------------------ */
/* Failures                                                                  */
/* ------------------------------------------------------------------------ */

static int fail_input(const subtraction *s, sub_input which, char *error,
                      size_t error_len)
{
    const char *why = h5reader_error(s->input[which]);

    snprintf(error, error_len, "%s: %s", s->path[which],
             why ? why : "unable to read it");
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

/* NaN in any input carries through every rule, so an unmeasured position and the columns
 * past a reference stay NaN. */

static void add_f32(const float *const *in, size_t n_in, float *out, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        float total = in[0][i];

        for (size_t k = 1; k < n_in; k++) {
            total += in[k][i];
        }

        out[i] = total;
    }
}

static void subtract_f32(const float *treated, const float *untreated, float *out,
                         size_t n)
{
    for (size_t i = 0; i < n; i++) {
        out[i] = treated[i] - untreated[i];
    }
}

/* Raises every negative value to zero, leaving NaN as it is. */
static void clip_f32(float *out, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        if (out[i] < 0.0F) {
            out[i] = 0.0F;
        }
    }
}

/* Each square is rounded before the two are added, which a single expression would let the
 * compiler fuse into a multiply-add. Fusing keeps one square at a wider precision than the
 * other, and the quadrature then depends on which input came first. */
static void propagate_f32(const float *treated, const float *untreated, float *out,
                          size_t n)
{
    for (size_t i = 0; i < n; i++) {
        float t = treated[i] * treated[i];
        float u = untreated[i] * untreated[i];

        out[i] = sqrtf(t + u);
    }
}

static void ratio_f32(const float *treated, const float *untreated,
                      const float *denatured, float *out, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        out[i] = denatured[i] > 0.0F
               ? (treated[i] - untreated[i]) / denatured[i]
               : (float)NAN;
    }
}

/* The error of the ratio, treating the three runs as independent.
 *
 * Not the relative-error form: that divides by the difference, which is zero at an
 * unreactive position, and squaring the control into a denominator underflows a float long
 * before the control does. */
static void ratio_error_f32(const float *const *rate, const float *const *error,
                            float *out, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        float d = rate[SUB_DENATURED][i];
        float r = (rate[SUB_TREATED][i] - rate[SUB_UNTREATED][i]) / d;
        float t = error[SUB_TREATED][i] * error[SUB_TREATED][i];
        float u = error[SUB_UNTREATED][i] * error[SUB_UNTREATED][i];
        float c = r * error[SUB_DENATURED][i];

        out[i] = d > 0.0F ? sqrtf(t + u + c * c) / d : (float)NAN;
    }
}

static void add_u64(const uint64_t *const *in, size_t n_in, uint64_t *out, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        uint64_t total = in[0][i];

        for (size_t k = 1; k < n_in; k++) {
            total += in[k][i];
        }

        out[i] = total;
    }
}

static void gather_f32(const subtraction *s, out_field_id id, const float **v)
{
    for (size_t i = 0; i < s->n_inputs; i++) {
        v[i] = s->rows[i].field[id];
    }
}

static int combine_f32(const subtraction *s, sub_rule how, const void *const *in,
                       float *out, size_t n)
{
    const float *v[SUB_N_INPUTS];
    const float *rate[SUB_N_INPUTS];

    for (size_t k = 0; k < s->n_inputs; k++) {
        v[k] = in[k];
    }

    switch (how) {
        case SUB_SUM:
            add_f32(v, s->n_inputs, out, n);
            return 0;
        case SUB_DIFFERENCE:
            subtract_f32(v[SUB_TREATED], v[SUB_UNTREATED], out, n);
            if (s->cfg->clip) {
                clip_f32(out, n);
            }
            return 0;
        case SUB_QUADRATURE:
            propagate_f32(v[SUB_TREATED], v[SUB_UNTREATED], out, n);
            return 0;
        case SUB_RATIO:
            ratio_f32(v[SUB_TREATED], v[SUB_UNTREATED], v[SUB_DENATURED], out, n);
            if (s->cfg->clip) {
                clip_f32(out, n);
            }
            return 0;
        case SUB_RATIO_ERROR:
            gather_f32(s, OUT_REACTIVITY, rate);
            ratio_error_f32(rate, v, out, n);
            return 0;
        case SUB_N_RULES:
            break;
    }

    return -1;
}

static int combine_u64(const subtraction *s, sub_rule how, const void *const *in,
                       uint64_t *out, size_t n)
{
    const uint64_t *v[SUB_N_INPUTS];

    if (how != SUB_SUM) {
        return -1;
    }

    for (size_t k = 0; k < s->n_inputs; k++) {
        v[k] = in[k];
    }

    add_u64(v, s->n_inputs, out, n);
    return 0;
}

/* in holds the field's own values, one from each input; a rule reading any other field
 * takes it from the subtraction. */
static int combine_values(const subtraction *s, out_field_id id,
                          const void *const *in, void *out, size_t n)
{
    sub_rule how = RULES[id][s->mode];

    switch (OUT_FIELDS[id].stored) {
        case OUT_F32:      return combine_f32(s, how, in, out, n);
        case OUT_U64:      return combine_u64(s, how, in, out, n);
        case OUT_N_STORED: break;
    }

    return -1;
}

/* ------------------------------------------------------------------------ */
/* Rows                                                                      */
/* ------------------------------------------------------------------------ */

static int read_reference(subtraction *s, int32_t tid, char *error, size_t error_len)
{
    for (size_t i = 0; i < s->n_inputs; i++) {
        for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
            if (!OUT_FIELDS[id].per_ref) {
                continue;
            }

            if (h5reader_field(s->input[i], id, tid, s->rows[i].field[id]) < 0) {
                return fail_input(s, (sub_input)i, error, error_len);
            }
        }
    }

    return 0;
}

static int combine_field(const subtraction *s, out_field_id id, size_t width)
{
    const void *in[SUB_N_INPUTS];

    for (size_t i = 0; i < s->n_inputs; i++) {
        in[i] = s->rows[i].field[id];
    }

    return combine_values(s, id, in, s->result, width);
}

static int write_field(subtraction *s, out_field_id id, int32_t tid, char *error,
                       size_t error_len)
{
    size_t width = out_values(id, s->ref_cap, s->ref_cap);

    if (combine_field(s, id, width) < 0) {
        snprintf(error, error_len, "%s: no rule combines a field of this type",
                 OUT_FIELDS[id].name);
        return -1;
    }

    /* Written whole rather than to the reference's own length: the columns past a
     * reference are NaN in every input and stay NaN through the arithmetic, so the writer
     * need not mark them. */
    if (h5writer_row(s->out, id, tid, s->result) < 0) {
        return fail_output(s, error, error_len);
    }

    return 0;
}

static int subtract_reference(subtraction *s, int32_t tid, char *error,
                              size_t error_len)
{
    if (read_reference(s, tid, error, error_len) < 0) {
        return -1;
    }

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (!OUT_FIELDS[id].per_ref) {
            continue;
        }

        if (write_field(s, id, tid, error, error_len) < 0) {
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

/* Combined before the row sets exist, so a field with no row may only follow a rule reading
 * the values passed to it. */
static int read_totals(subtraction *s, char *error, size_t error_len)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        uint64_t    held[SUB_N_INPUTS];
        const void *in[SUB_N_INPUTS];
        uint64_t    combined;

        if (OUT_FIELDS[id].per_ref) {
            continue;
        }

        for (size_t i = 0; i < s->n_inputs; i++) {
            size_t total;

            if (h5reader_total(s->input[i], id, &total) < 0) {
                return fail_input(s, (sub_input)i, error, error_len);
            }

            held[i] = total;
            in[i]   = &held[i];
        }

        if (combine_values(s, id, in, &combined, 1) < 0) {
            snprintf(error, error_len, "%s: no rule combines a field of this type",
                     OUT_FIELDS[id].name);
            return -1;
        }

        s->totals[id] = (size_t)combined;
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

static int check_agreement(subtraction *s, char *error, size_t error_len)
{
    for (size_t i = 1; i < s->n_inputs; i++) {
        if (h5reader_refs(s->input[i]) != h5reader_refs(s->input[SUB_TREATED])) {
            snprintf(error, error_len,
                     "%s holds %d references and %s holds %d",
                     s->path[SUB_TREATED], h5reader_refs(s->input[SUB_TREATED]),
                     s->path[i], h5reader_refs(s->input[i]));
            return -1;
        }

        if (h5reader_capacity(s->input[i]) != h5reader_capacity(s->input[SUB_TREATED])) {
            snprintf(error, error_len,
                     "%s is %zu bases wide and %s is %zu; the two were counted "
                     "against different references",
                     s->path[SUB_TREATED], h5reader_capacity(s->input[SUB_TREATED]),
                     s->path[i], h5reader_capacity(s->input[i]));
            return -1;
        }
    }

    s->n_refs  = h5reader_refs(s->input[SUB_TREATED]);
    s->ref_cap = h5reader_capacity(s->input[SUB_TREATED]);

    return 0;
}

static void take_inputs(subtraction *s)
{
    s->path[SUB_TREATED]   = s->cfg->treated_path;
    s->path[SUB_UNTREATED] = s->cfg->untreated_path;
    s->n_inputs            = 2;
    s->mode                = SUB_UNCONTROLLED;

    if (s->cfg->denatured_path) {
        s->path[SUB_DENATURED] = s->cfg->denatured_path;
        s->n_inputs            = 3;
        s->mode                = SUB_CONTROLLED;
    }
}

static int open_inputs(subtraction *s, char *error, size_t error_len)
{
    take_inputs(s);

    for (size_t i = 0; i < s->n_inputs; i++) {
        s->input[i] = h5reader_open(s->path[i]);

        if (!s->input[i]) {
            snprintf(error, error_len, "out of memory");
            return -1;
        }
    }

    for (size_t i = 0; i < s->n_inputs; i++) {
        if (h5reader_error(s->input[i])) {
            return fail_input(s, (sub_input)i, error, error_len);
        }
    }

    if (check_agreement(s, error, error_len) < 0) {
        return -1;
    }

    return read_totals(s, error, error_len);
}

static int build_buffers(subtraction *s, char *error, size_t error_len)
{
    for (size_t i = 0; i < s->n_inputs; i++) {
        for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
            if (!OUT_FIELDS[id].per_ref) {
                continue;
            }

            s->rows[i].field[id] = calloc(out_values(id, s->ref_cap, s->ref_cap),
                                          out_stored_bytes(id));

            if (!s->rows[i].field[id]) {
                snprintf(error, error_len, "out of memory");
                return -1;
            }
        }
    }

    s->result = calloc(out_widest(s->ref_cap), out_widest_bytes());

    if (!s->result) {
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

    for (size_t i = 0; i < SUB_N_INPUTS; i++) {
        h5reader_close(s->input[i]);

        for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
            free(s->rows[i].field[id]);
        }
    }

    free(s->result);
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

    /* Every input is checked before the output is created, so a run that refuses its
     * inputs leaves the file at that path alone. */
    if (open_inputs(&s, error, error_len) == 0 &&
        build_buffers(&s, error, error_len) == 0 &&
        open_output(&s, may_replace, error, error_len) == 0 &&
        subtract_references(&s, error, error_len) == 0) {
        status = write_totals(&s, error, error_len);
    }

    subtraction_teardown(&s);
    return status;
}

/* ------------------------------------------------------------------------ */
/* The rules, described                                                      */
/* ------------------------------------------------------------------------ */

void subtract_dump_rules(FILE *out)
{
    fprintf(out, "{\n  \"version\": \"%s\",\n  \"fields\": [\n", CMUTS_VERSION);

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        const sub_rule_description *alone = &DESCRIPTIONS[RULES[id][SUB_UNCONTROLLED]];
        const sub_rule_description *controlled = &DESCRIPTIONS[RULES[id][SUB_CONTROLLED]];

        fprintf(out,
                "    {\n"
                "      \"name\": \"%s\",\n"
                "      \"uncontrolled\": { \"rule\": \"%s\", \"detail\": \"%s\" },\n"
                "      \"controlled\": { \"rule\": \"%s\", \"detail\": \"%s\" }\n"
                "    }%s\n",
                OUT_FIELDS[id].name,
                alone->name, alone->detail,
                controlled->name, controlled->detail,
                id + 1 < OUT_N_FIELDS ? "," : "");
    }

    fprintf(out, "  ]\n}\n");
}
