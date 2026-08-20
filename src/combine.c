/* combine.c -- combining several outputs into one, field by field.
 *
 * Every input is read a row at a time and the result written the same way, so memory is
 * bounded by the longest reference and not by the size of the files.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "combine.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "format.h"
#include "h5reader.h"
#include "h5writer.h"
#include "progress.h"

struct combine_rows {
    void **value;  /* n_inputs * FMT_N_FIELDS buffers, indexed by input then field */
    size_t n_inputs;
};

typedef struct {
    const combine_spec *spec;

    h5reader **input;
    h5writer  *out;
    fmt_request requests[FMT_N_FIELDS];  /* derived from the spec's manifest */
    size_t      n_requests;

    combine_rows rows;
    void        *result;

    int32_t n_refs;
    size_t  ref_cap;
    bool    writes[FMT_N_FIELDS];   /* what the spec says this run leaves behind */
} combination;

/* ------------------------------------------------------------------------ */
/* Rows                                                                      */
/* ------------------------------------------------------------------------ */

static void *row_at(const combine_rows *rows, size_t input, fmt_field_id id)
{
    return rows->value[(input * FMT_N_FIELDS) + id];
}

const void *combine_row(const combine_rows *rows, size_t input, fmt_field_id id)
{
    return row_at(rows, input, id);
}

/* Gives the values one row of a field holds: one for a field with no row, and none for a
 * field this program does not write, so no buffer is taken for it. */
static size_t row_values(const combination *c, fmt_field_id id, size_t ref_cap)
{
    if (!fmt_wanted(id, c->writes)) {
        return 0;
    }

    return FMT_FIELDS[id].per_ref ? fmt_values(id, ref_cap, ref_cap) : 1;
}

/* ------------------------------------------------------------------------ */
/* Failures                                                                  */
/* ------------------------------------------------------------------------ */

static int fail_input(const combination *c, size_t which, char *error, size_t error_len)
{
    return h5reader_fail(c->input[which], c->spec->inputs[which], error, error_len);
}

static int fail_output(const combination *c, char *error, size_t error_len)
{
    return h5writer_fail(c->out, c->spec->output, error, error_len);
}

static int fail_rule(fmt_field_id id, int status, char *error, size_t error_len)
{
    const char *why = status == COMBINE_MISMATCH
                    ? "the inputs disagree on it, so they were not made against one FASTA"
                    : "no rule combines a field of this type";

    snprintf(error, error_len, "%s: %s", FMT_FIELDS[id].name, why);
    return -1;
}

/* ------------------------------------------------------------------------ */
/* Summing                                                                   */
/* ------------------------------------------------------------------------ */

/* Accumulating in a local and storing once rounds only at the end, so the sum does not
 * depend on the width the platform evaluates floats at. */

static void sum_f32(const combine_rows *rows, fmt_field_id id, float *out, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        float total = ((const float *)row_at(rows, 0, id))[i];

        for (size_t k = 1; k < rows->n_inputs; k++) {
            total += ((const float *)row_at(rows, k, id))[i];
        }

        out[i] = total;
    }
}

static void sum_u64(const combine_rows *rows, fmt_field_id id, uint64_t *out, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        uint64_t total = ((const uint64_t *)row_at(rows, 0, id))[i];

        for (size_t k = 1; k < rows->n_inputs; k++) {
            total += ((const uint64_t *)row_at(rows, k, id))[i];
        }

        out[i] = total;
    }
}

int combine_sum(const combine_rows *rows, fmt_field_id id, void *out, size_t n)
{
    switch (FMT_FIELDS[id].stored) {
        case FMT_F32:
            sum_f32(rows, id, out, n);
            return 0;
        case FMT_U64:
            sum_u64(rows, id, out, n);
            return 0;
        case FMT_I8:
        case FMT_N_STORED:
            break;
    }

    return COMBINE_NO_RULE;
}

int combine_same(const combine_rows *rows, fmt_field_id id, void *out, size_t n)
{
    size_t      bytes = n * fmt_stored_bytes(id);
    const void *first = combine_row(rows, 0, id);

    for (size_t i = 1; i < rows->n_inputs; i++) {
        if (memcmp(first, combine_row(rows, i, id), bytes) != 0) {
            return COMBINE_MISMATCH;
        }
    }

    memcpy(out, first, bytes);
    return COMBINE_OK;
}

/* ------------------------------------------------------------------------ */
/* References                                                                */
/* ------------------------------------------------------------------------ */

static int read_reference(combination *c, int32_t tid, char *error, size_t error_len)
{
    for (size_t i = 0; i < c->spec->n_inputs; i++) {
        for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
            if (!FMT_FIELDS[id].per_ref || !fmt_wanted(id, c->writes)) {
                continue;
            }

            if (h5reader_field(c->input[i], id, tid, row_at(&c->rows, i, id)) < 0) {
                return fail_input(c, i, error, error_len);
            }
        }
    }

    return 0;
}

static int combine_field(const combination *c, fmt_field_id id, size_t n, char *error,
                         size_t error_len)
{
    int status = c->spec->field(&c->rows, id, c->result, n, c->spec->ctx);

    if (status < 0) {
        return fail_rule(id, status, error, error_len);
    }

    return 0;
}

static int write_reference(combination *c, int32_t tid, char *error, size_t error_len)
{
    for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
        if (!FMT_FIELDS[id].per_ref || !fmt_wanted(id, c->writes)) {
            continue;
        }

        if (combine_field(c, id, row_values(c, id, c->ref_cap), error, error_len) < 0) {
            return -1;
        }

        /* Written whole and not to the reference's own length: the columns past a
         * reference are NaN in every input and stay NaN through the arithmetic, so the
         * writer does not mark them. */
        if (h5writer_row(c->out, id, tid, c->result) < 0) {
            return fail_output(c, error, error_len);
        }
    }

    return 0;
}

static int combine_references(combination *c, char *error, size_t error_len)
{
    progress *bar    = progress_start((uint64_t)c->n_refs);
    int       status = 0;

    for (int32_t tid = 0; tid < c->n_refs && status == 0; tid++) {
        if (read_reference(c, tid, error, error_len) < 0 ||
            write_reference(c, tid, error, error_len) < 0) {
            status = -1;
        }

        progress_follow(bar, (uint64_t)tid + 1);
    }

    progress_finish(bar);
    return status;
}

/* ------------------------------------------------------------------------ */
/* Totals                                                                    */
/* ------------------------------------------------------------------------ */

/* A field with no row holds one value per input, read once into the row set beside the
 * rows. h5reader gives it as a count, so it is stored as a whole unsigned. */

static int read_totals(combination *c, char *error, size_t error_len)
{
    for (size_t i = 0; i < c->spec->n_inputs; i++) {
        for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
            size_t total;

            if (FMT_FIELDS[id].per_ref || !fmt_wanted(id, c->writes)) {
                continue;
            }

            if (h5reader_total(c->input[i], id, &total) < 0) {
                return fail_input(c, i, error, error_len);
            }

            *(uint64_t *)row_at(&c->rows, i, id) = total;
        }
    }

    return 0;
}

static int write_totals(combination *c, char *error, size_t error_len)
{
    for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
        if (FMT_FIELDS[id].per_ref || !fmt_wanted(id, c->writes)) {
            continue;
        }

        if (combine_field(c, id, 1, error, error_len) < 0) {
            return -1;
        }

        if (h5writer_total(c->out, id, (size_t) * (uint64_t *)c->result) < 0) {
            return fail_output(c, error, error_len);
        }
    }

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Assembly                                                                  */
/* ------------------------------------------------------------------------ */

static int check_agreement(combination *c, char *error, size_t error_len)
{
    for (size_t i = 1; i < c->spec->n_inputs; i++) {
        if (h5reader_refs(c->input[i]) != h5reader_refs(c->input[0])) {
            snprintf(error, error_len, "%s holds %d references and %s holds %d",
                     c->spec->inputs[0], h5reader_refs(c->input[0]),
                     c->spec->inputs[i], h5reader_refs(c->input[i]));
            return -1;
        }

        if (h5reader_capacity(c->input[i]) != h5reader_capacity(c->input[0])) {
            snprintf(error, error_len,
                     "%s is %zu bases wide and %s is %zu; the two were counted "
                     "against different references",
                     c->spec->inputs[0], h5reader_capacity(c->input[0]),
                     c->spec->inputs[i], h5reader_capacity(c->input[i]));
            return -1;
        }
    }

    c->n_refs  = h5reader_refs(c->input[0]);
    c->ref_cap = h5reader_capacity(c->input[0]);

    return 0;
}

/* Clears from writes every field that depends on one some input does not carry. What
 * remains is read, allocated for and written alike. */
static void drop_absent_fields(combination *c)
{
    for (size_t i = 0; i < c->spec->writes->n_fields; i++) {
        const fmt_written *field = &c->spec->writes->fields[i];

        for (const fmt_field_id *dep = field->depends;
             dep && *dep != FMT_N_FIELDS; dep++) {
            for (size_t k = 0; k < c->spec->n_inputs; k++) {
                if (!h5reader_holds(c->input[k], *dep)) {
                    c->writes[field->id] = false;
                }
            }
        }
    }
}

static int open_inputs(combination *c, char *error, size_t error_len)
{
    c->input = calloc(c->spec->n_inputs, sizeof *c->input);

    if (!c->input) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    c->n_requests = fmt_requests_of(c->spec->writes, c->requests);

    for (size_t i = 0; i < c->spec->n_inputs; i++) {
        c->input[i] = h5reader_open(c->spec->inputs[i], c->requests, c->n_requests);

        if (!c->input[i]) {
            snprintf(error, error_len, "out of memory");
            return -1;
        }
    }

    for (size_t i = 0; i < c->spec->n_inputs; i++) {
        if (h5reader_error(c->input[i])) {
            return fail_input(c, i, error, error_len);
        }
    }

    for (size_t i = 0; i < c->spec->n_inputs; i++) {
        const char *ignored = h5reader_ignored(c->input[i]);

        if (ignored) {
            fprintf(stderr, "%s: ignoring %s\n", c->spec->inputs[i], ignored);
        }
    }

    if (check_agreement(c, error, error_len) < 0) {
        return -1;
    }

    drop_absent_fields(c);
    return 0;
}

static int build_rows(combination *c, char *error, size_t error_len)
{
    c->rows.n_inputs = c->spec->n_inputs;
    c->rows.value    = calloc(c->spec->n_inputs * FMT_N_FIELDS, sizeof *c->rows.value);

    if (!c->rows.value) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    for (size_t i = 0; i < c->spec->n_inputs; i++) {
        for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
            void **slot = &c->rows.value[(i * FMT_N_FIELDS) + id];

            size_t values = row_values(c, id, c->ref_cap);

            if (values == 0) {
                continue;
            }

            *slot = calloc(values, fmt_stored_bytes(id));

            if (!*slot) {
                snprintf(error, error_len, "out of memory");
                return -1;
            }
        }
    }

    c->result = calloc(fmt_widest(c->ref_cap, c->writes), fmt_widest_bytes());

    if (!c->result) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    return 0;
}

static int open_output(combination *c, bool may_replace, char *error, size_t error_len)
{
    c->out = h5writer_create(c->spec->output, c->spec->program, c->n_refs, c->ref_cap,
                             may_replace, c->writes, false);
    if (!c->out) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    return h5writer_error(c->out) ? fail_output(c, error, error_len) : 0;
}

static void combination_teardown(combination *c)
{
    h5writer_close(c->out);

    if (c->input) {
        for (size_t i = 0; i < c->spec->n_inputs; i++) {
            h5reader_close(c->input[i]);
        }
    }

    if (c->rows.value) {
        for (size_t i = 0; i < c->spec->n_inputs * FMT_N_FIELDS; i++) {
            free(c->rows.value[i]);
        }
    }

    free(c->rows.value);
    free(c->input);
    free(c->result);
}

int combine_run(const combine_spec *spec, char *error, size_t error_len)
{
    combination c           = { .spec = spec };

    fmt_selection(spec->writes, c.writes);
    bool        may_replace = false;
    int         status      = -1;

    if (h5writer_may_replace(spec->output, spec->overwrite, &may_replace, error,
                             error_len) < 0) {
        return -1;
    }

    /* Every input is checked before the output is created, so a run that fails on its
     * inputs does not modify the file at that path. */
    if (open_inputs(&c, error, error_len) == 0 &&
        build_rows(&c, error, error_len) == 0 &&
        read_totals(&c, error, error_len) == 0 &&
        open_output(&c, may_replace, error, error_len) == 0 &&
        combine_references(&c, error, error_len) == 0) {
        status = write_totals(&c, error, error_len);
    }

    combination_teardown(&c);
    return status;
}
