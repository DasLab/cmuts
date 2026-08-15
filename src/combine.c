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

#include "h5reader.h"
#include "h5writer.h"
#include "output.h"

struct combine_rows {
    void **value;  /* n_inputs * OUT_N_FIELDS buffers, indexed by input then field */
    size_t n_inputs;
};

typedef struct {
    const combine_spec *spec;

    h5reader **input;
    h5writer  *out;

    combine_rows rows;
    void        *result;

    int32_t n_refs;
    size_t  ref_cap;
} combination;

/* ------------------------------------------------------------------------ */
/* Rows                                                                      */
/* ------------------------------------------------------------------------ */

static void *row_at(const combine_rows *rows, size_t input, out_field_id id)
{
    return rows->value[(input * OUT_N_FIELDS) + id];
}

const void *combine_row(const combine_rows *rows, size_t input, out_field_id id)
{
    return row_at(rows, input, id);
}

/* Gives the values one field occupies in a row set. A field with no row holds one
 * value. */
static size_t row_values(out_field_id id, size_t ref_cap)
{
    return OUT_FIELDS[id].per_ref ? out_values(id, ref_cap, ref_cap) : 1;
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

static int fail_rule(out_field_id id, char *error, size_t error_len)
{
    snprintf(error, error_len, "%s: no rule combines a field of this type",
             OUT_FIELDS[id].name);
    return -1;
}

/* ------------------------------------------------------------------------ */
/* Summing                                                                   */
/* ------------------------------------------------------------------------ */

/* Accumulating in a local and storing once rounds only at the end, whatever width the
 * platform evaluates floats at. */

static void sum_f32(const combine_rows *rows, out_field_id id, float *out, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        float total = ((const float *)row_at(rows, 0, id))[i];

        for (size_t k = 1; k < rows->n_inputs; k++) {
            total += ((const float *)row_at(rows, k, id))[i];
        }

        out[i] = total;
    }
}

static void sum_u64(const combine_rows *rows, out_field_id id, uint64_t *out, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        uint64_t total = ((const uint64_t *)row_at(rows, 0, id))[i];

        for (size_t k = 1; k < rows->n_inputs; k++) {
            total += ((const uint64_t *)row_at(rows, k, id))[i];
        }

        out[i] = total;
    }
}

int combine_sum(const combine_rows *rows, out_field_id id, void *out, size_t n)
{
    switch (OUT_FIELDS[id].stored) {
        case OUT_F32:
            sum_f32(rows, id, out, n);
            return 0;
        case OUT_U64:
            sum_u64(rows, id, out, n);
            return 0;
        case OUT_N_STORED:
            break;
    }

    return -1;
}

/* ------------------------------------------------------------------------ */
/* References                                                                */
/* ------------------------------------------------------------------------ */

static int read_reference(combination *c, int32_t tid, char *error, size_t error_len)
{
    for (size_t i = 0; i < c->spec->n_inputs; i++) {
        for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
            if (!OUT_FIELDS[id].per_ref) {
                continue;
            }

            if (h5reader_field(c->input[i], id, tid, row_at(&c->rows, i, id)) < 0) {
                return fail_input(c, i, error, error_len);
            }
        }
    }

    return 0;
}

static int combine_field(const combination *c, out_field_id id, size_t n, char *error,
                         size_t error_len)
{
    if (c->spec->field(&c->rows, id, c->result, n, c->spec->ctx) < 0) {
        return fail_rule(id, error, error_len);
    }

    return 0;
}

static int write_reference(combination *c, int32_t tid, char *error, size_t error_len)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (!OUT_FIELDS[id].per_ref) {
            continue;
        }

        if (combine_field(c, id, row_values(id, c->ref_cap), error, error_len) < 0) {
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
    for (int32_t tid = 0; tid < c->n_refs; tid++) {
        if (read_reference(c, tid, error, error_len) < 0 ||
            write_reference(c, tid, error, error_len) < 0) {
            return -1;
        }
    }

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Totals                                                                    */
/* ------------------------------------------------------------------------ */

/* A field with no row holds one value per input, read once into the row set beside the
 * rows. h5reader gives it as a count, so it is stored as a whole unsigned. */

static int read_totals(combination *c, char *error, size_t error_len)
{
    for (size_t i = 0; i < c->spec->n_inputs; i++) {
        for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
            size_t total;

            if (OUT_FIELDS[id].per_ref) {
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
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        if (OUT_FIELDS[id].per_ref) {
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

static int open_inputs(combination *c, char *error, size_t error_len)
{
    c->input = calloc(c->spec->n_inputs, sizeof *c->input);

    if (!c->input) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    for (size_t i = 0; i < c->spec->n_inputs; i++) {
        c->input[i] = h5reader_open(c->spec->inputs[i]);

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

    return check_agreement(c, error, error_len);
}

static int build_rows(combination *c, char *error, size_t error_len)
{
    c->rows.n_inputs = c->spec->n_inputs;
    c->rows.value    = calloc(c->spec->n_inputs * OUT_N_FIELDS, sizeof *c->rows.value);

    if (!c->rows.value) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    for (size_t i = 0; i < c->spec->n_inputs; i++) {
        for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
            void **slot = &c->rows.value[(i * OUT_N_FIELDS) + id];

            *slot = calloc(row_values(id, c->ref_cap), out_stored_bytes(id));

            if (!*slot) {
                snprintf(error, error_len, "out of memory");
                return -1;
            }
        }
    }

    c->result = calloc(out_widest(c->ref_cap), out_widest_bytes());

    if (!c->result) {
        snprintf(error, error_len, "out of memory");
        return -1;
    }

    return 0;
}

static int open_output(combination *c, bool may_replace, char *error, size_t error_len)
{
    c->out = h5writer_create(c->spec->output, c->spec->program, c->n_refs, c->ref_cap,
                             may_replace);
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
        for (size_t i = 0; i < c->spec->n_inputs * OUT_N_FIELDS; i++) {
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
