/* csv.c -- a cmuts output written as comma separated values.
 *
 * The sequence dataset supplies the base at each position and the length of each
 * reference, so the input alone determines the rows. A FASTA supplies the names, and is
 * checked against the input. Without one, the rows are numbered.
 *
 * Each requested field the input holds becomes one column, under the name of its dataset
 * in format.h. One reference is read at a time, so memory is bounded by the longest
 * reference and not by the size of the input.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "csv.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fasta.h"
#include "format.h"
#include "h5reader.h"
#include "nuc.h"

/* Characters that force a value to be quoted: the separator, the quote itself, and the
 * two line endings. */
#define QUOTED_BY ",\"\r\n"

/* The columns identifying a row, written before those read from the input. */
#define KEY_COLUMNS "reference,position,base"

/* Room for a row number, which is the reference name when no FASTA is given. */
#define NUMBER_MAX 32

/* ------------------------------------------------------------------------ */
/* The run                                                                   */
/* ------------------------------------------------------------------------ */

/* The state of one run. columns and values are parallel: values[i] holds one reference's
 * row of columns[i]. The sequence is read into bases instead, because it supplies the base
 * key column and not a column of its own. */
typedef struct {
    const csv_config *cfg;
    FILE             *out;
    h5reader         *reader;
    cm_fasta_reader  *fasta;   /* NULL where no FASTA was given */
    fmt_field_id      columns[FMT_N_FIELDS];
    float            *values[FMT_N_FIELDS];
    int8_t           *bases;
    size_t            n_columns;
    size_t            cap;
} context;

/* Selects one column for each requested field the reader opened, in the order of the
 * requests. The sequence is excluded, because it supplies the base key column. */
static void columns_of(context *ctx, const fmt_request *requests, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        if (requests[i].id != FMT_SEQUENCE && h5reader_holds(ctx->reader, requests[i].id)) {
            ctx->columns[ctx->n_columns++] = requests[i].id;
        }
    }
}

static int fail_memory(char *error, size_t error_len)
{
    snprintf(error, error_len, "out of memory");
    return -1;
}

static int allocate(context *ctx, size_t cap, char *error, size_t error_len)
{
    ctx->cap   = cap;
    ctx->bases = malloc(fmt_values(FMT_SEQUENCE, cap, cap) * sizeof *ctx->bases);

    if (!ctx->bases) {
        return fail_memory(error, error_len);
    }

    for (size_t i = 0; i < ctx->n_columns; i++) {
        size_t values = fmt_values(ctx->columns[i], cap, cap);

        ctx->values[i] = malloc(values * sizeof *ctx->values[i]);

        if (!ctx->values[i]) {
            return fail_memory(error, error_len);
        }
    }

    return 0;
}

static void release(context *ctx)
{
    for (size_t i = 0; i < ctx->n_columns; i++) {
        free(ctx->values[i]);
    }

    free(ctx->bases);
}

/* ------------------------------------------------------------------------ */
/* The table                                                                 */
/* ------------------------------------------------------------------------ */

/* Prints one text value. A value holding any character of QUOTED_BY is printed in quotes,
 * with each of its own quotes doubled. */
static void print_field(FILE *out, const char *text)
{
    if (!strpbrk(text, QUOTED_BY)) {
        fputs(text, out);
        return;
    }

    fputc('"', out);

    for (const char *at = text; *at; at++) {
        if (*at == '"') {
            fputc('"', out);
        }

        fputc(*at, out);
    }

    fputc('"', out);
}

/* Prints one value of a column. A value the run did not measure is not a number, and is
 * printed as an empty field. */
static void print_value(FILE *out, float value)
{
    if (isfinite(value)) {
        fprintf(out, "%.6g", (double)value);
    }
}

/* Prints the column names: KEY_COLUMNS, then one for each field read from the input. */
static void print_header(const context *ctx)
{
    fputs(KEY_COLUMNS, ctx->out);

    for (size_t i = 0; i < ctx->n_columns; i++) {
        fprintf(ctx->out, ",%s", FMT_FIELDS[ctx->columns[i]].name);
    }

    fputc('\n', ctx->out);
}

/* Prints one row: the reference, the position counted from one, the base there, and each
 * column's value there. */
static void print_position(const context *ctx, const char *name, size_t at)
{
    print_field(ctx->out, name);
    fprintf(ctx->out, ",%zu,%c", at + 1, nuc_char(ctx->bases[at]));

    for (size_t i = 0; i < ctx->n_columns; i++) {
        fputc(',', ctx->out);
        print_value(ctx->out, ctx->values[i][at]);
    }

    fputc('\n', ctx->out);
}

/* ------------------------------------------------------------------------ */
/* One reference                                                             */
/* ------------------------------------------------------------------------ */

static int fail_read(const context *ctx, char *error, size_t error_len)
{
    snprintf(error, error_len, "%s: %s", ctx->cfg->input_path,
             h5reader_error(ctx->reader));

    return -1;
}

/* Reads one reference's sequence and its row of every column. */
static int read_row(context *ctx, int32_t tid, char *error, size_t error_len)
{
    for (size_t i = 0; i < ctx->n_columns; i++) {
        if (h5reader_field(ctx->reader, ctx->columns[i], tid, ctx->values[i]) < 0) {
            return fail_read(ctx, error, error_len);
        }
    }

    if (h5reader_field(ctx->reader, FMT_SEQUENCE, tid, ctx->bases) < 0) {
        return fail_read(ctx, error, error_len);
    }

    return 0;
}

/* Returns the length of one reference, which is the number of values before the first
 * padding value. */
static size_t length_of(const int8_t *bases, size_t cap)
{
    size_t len = 0;

    while (len < cap && bases[len] != NUC_TOKEN_NONE) {
        len++;
    }

    return len;
}

/* Fills error where the row is not a run of tokens followed by padding, which is the only
 * shape a run writes. Every value before len is not padding, since len ends there, so a
 * value there need only be checked for being a token. */
static int check_row(const context *ctx, int32_t tid, size_t len, char *error,
                     size_t error_len)
{
    for (size_t at = 0; at < ctx->cap; at++) {
        bool valid = at < len ? nuc_is_token(ctx->bases[at])
                              : ctx->bases[at] == NUC_TOKEN_NONE;

        if (!valid) {
            snprintf(error, error_len,
                     "%s: row %d is not a sequence; position %zu holds %d",
                     ctx->cfg->input_path, tid + 1, at + 1, ctx->bases[at]);
            return -1;
        }
    }

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Naming the rows                                                           */
/* ------------------------------------------------------------------------ */

/* Fills error with the FASTA reader's failure, and returns whether it had one. */
static bool fail_fasta(const context *ctx, char *error, size_t error_len)
{
    if (!cm_fasta_error(ctx->fasta)) {
        return false;
    }

    snprintf(error, error_len, "%s: %s", ctx->cfg->fasta_path,
             cm_fasta_error(ctx->fasta));

    return true;
}

/* Fills error where the FASTA cannot be read, or holds fewer records than the input holds
 * rows. */
static void fail_names(const context *ctx, char *error, size_t error_len)
{
    if (fail_fasta(ctx, error, error_len)) {
        return;
    }

    snprintf(error, error_len, "%s names fewer references than the %d rows of %s",
             ctx->cfg->fasta_path, h5reader_refs(ctx->reader), ctx->cfg->input_path);
}

/* Returns the name of one row: the name of the next FASTA record, or the row number
 * written into number. Returns NULL and fills error where the FASTA does not describe this
 * row. A record's name borrows the reader's memory, so the caller must use it before
 * asking for the next row.
 *
 * A record of another length describes another molecule. Such a FASTA is refused, because
 * otherwise every row after the disagreement would carry the wrong name. */
static const char *name_of(context *ctx, int32_t tid, size_t len, char *number,
                           size_t number_len, char *error, size_t error_len)
{
    cm_fasta_record record;

    if (!ctx->fasta) {
        snprintf(number, number_len, "%d", tid + 1);
        return number;
    }

    if (cm_fasta_next(ctx->fasta, &record) != CM_ITER_OK) {
        fail_names(ctx, error, error_len);
        return NULL;
    }

    if (record.len != len) {
        snprintf(error, error_len,
                 "%s is %zu long and row %d of %s is %zu; the FASTA is not the one it was "
                 "counted against", record.name, record.len, tid + 1,
                 ctx->cfg->input_path, len);
        return NULL;
    }

    return record.name;
}

/* Fails where the FASTA holds more records than the input holds rows. */
static int check_names_end(context *ctx, char *error, size_t error_len)
{
    cm_fasta_record record;

    if (!ctx->fasta) {
        return 0;
    }

    if (cm_fasta_next(ctx->fasta, &record) == CM_ITER_OK) {
        snprintf(error, error_len, "%s names more references than the %d rows of %s",
                 ctx->cfg->fasta_path, h5reader_refs(ctx->reader), ctx->cfg->input_path);
        return -1;
    }

    return fail_fasta(ctx, error, error_len) ? -1 : 0;
}

/* ------------------------------------------------------------------------ */
/* The run                                                                   */
/* ------------------------------------------------------------------------ */

/* Writes every row of the input, in the order the input holds them. */
static int write_all(context *ctx, char *error, size_t error_len)
{
    char number[NUMBER_MAX];

    print_header(ctx);

    for (int32_t tid = 0; tid < h5reader_refs(ctx->reader); tid++) {
        const char *name;
        size_t      len;

        if (read_row(ctx, tid, error, error_len) < 0) {
            return -1;
        }

        len = length_of(ctx->bases, ctx->cap);

        if (check_row(ctx, tid, len, error, error_len) < 0) {
            return -1;
        }

        name = name_of(ctx, tid, len, number, sizeof number, error, error_len);

        if (!name) {
            return -1;
        }

        for (size_t at = 0; at < len; at++) {
            print_position(ctx, name, at);
        }
    }

    return check_names_end(ctx, error, error_len);
}

static int open_names(context *ctx, char *error, size_t error_len)
{
    const char *why = NULL;

    if (!ctx->cfg->fasta_path) {
        return 0;
    }

    ctx->fasta = cm_fasta_open(ctx->cfg->fasta_path, &why);

    if (!ctx->fasta) {
        snprintf(error, error_len, "%s: %s", ctx->cfg->fasta_path, why);
        return -1;
    }

    return 0;
}

int csv_run(const csv_config *cfg, const fmt_request *requests, size_t n_requests,
            FILE *out, char *error, size_t error_len)
{
    context ctx = { .cfg = cfg, .out = out };
    int     status;

    ctx.reader = h5reader_open(cfg->input_path, requests, n_requests);

    if (!ctx.reader) {
        return fail_memory(error, error_len);
    }

    if (h5reader_error(ctx.reader)) {
        status = h5reader_fail(ctx.reader, cfg->input_path, error, error_len);
        h5reader_close(ctx.reader);
        return status;
    }

    columns_of(&ctx, requests, n_requests);

    if (open_names(&ctx, error, error_len) < 0
        || allocate(&ctx, h5reader_capacity(ctx.reader), error, error_len) < 0) {
        status = -1;
    } else {
        status = write_all(&ctx, error, error_len);
    }

    release(&ctx);
    cm_fasta_close(ctx.fasta);
    h5reader_close(ctx.reader);

    return status;
}
