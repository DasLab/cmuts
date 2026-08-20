/* output.c -- the output field table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "output.h"

#include <math.h>

#include "version.h"

const out_field OUT_FIELDS[OUT_N_FIELDS] = {
    [OUT_COVERAGE] = {
        .name    = "coverage",
        .row     = shape_per_base,
        .per_ref = true,
        .stored  = OUT_F32,
        .fill    = 0.0,
    },
    [OUT_REACTIVITY] = {
        .name    = "reactivity",
        .row     = shape_per_base,
        .per_ref = true,
        .stored  = OUT_F32,
        .fill    = (double)NAN,
    },
    [OUT_ERROR] = {
        .name    = "error",
        .row     = shape_per_base,
        .per_ref = true,
        .stored  = OUT_F32,
        .fill    = (double)NAN,
    },
    [OUT_LENGTHS] = {
        .name    = "reads/lengths",
        .row     = shape_per_length,
        .per_ref = true,
        .stored  = OUT_U64,
        .fill    = 0.0,
    },
    [OUT_READS] = {
        .name    = "reads/counted",
        .row     = shape_none,
        .per_ref = true,
        .stored  = OUT_U64,
        .fill    = 0.0,
    },
    [OUT_REJECTED] = {
        .name    = "reads/rejected",
        .row     = shape_none,
        .per_ref = true,
        .stored  = OUT_U64,
        .fill    = 0.0,
    },
    [OUT_PAIRWISE_CORRELATION] = {
        .name     = "pairwise/correlation",
        .row      = shape_per_pair,
        .per_ref  = true,
        .stored   = OUT_F32,
        .fill     = (double)NAN,
    },
    [OUT_PAIRWISE_CONDITIONAL] = {
        .name     = "pairwise/conditional",
        .row      = shape_per_pair,
        .per_ref  = true,
        .stored   = OUT_F32,
        .fill     = (double)NAN,
    },
    [OUT_PAIRWISE_COVERAGE] = {
        .name     = "pairwise/coverage",
        .row      = shape_per_pair,
        .per_ref  = true,
        .stored   = OUT_F32,
        .fill     = 0.0,
    },
    [OUT_NORM] = {
        .name     = "norm",
        .row      = shape_none,
        .per_ref  = false,
        .stored   = OUT_F32,
        .fill     = (double)NAN,
    },
    [OUT_SEQUENCE] = {
        .name     = "sequence",
        .row      = shape_per_base,
        .per_ref  = true,
        .from_ref = true,
        .stored   = OUT_I8,
        .fill     = -1.0,
    },
    [OUT_UNMAPPED] = {
        .name    = "reads/unmapped",
        .row     = shape_none,
        .per_ref = false,
        .stored  = OUT_U64,
        .fill    = 0.0,
    },
};

const out_attribute OUT_ATTRIBUTES[OUT_N_ATTRS] = {
    [OUT_ATTR_PROGRAM] = {
        .name   = "program",
        .detail = "The name of the program that produced this file.",
    },
    [OUT_ATTR_VERSION] = {
        .name   = "version",
        .detail = "The version of cmuts that produced this file.",
    },
};

size_t out_values(out_field_id id, size_t len, size_t cap)
{
    return shape_values(OUT_FIELDS[id].row, len, cap);
}

void out_selection(const out_manifest *manifest, bool *wanted)
{
    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        wanted[id] = false;
    }

    for (size_t i = 0; i < manifest->n_fields; i++) {
        wanted[manifest->fields[i].id] = true;
    }
}

out_origin out_origin_of(const out_manifest *manifest, out_field_id id)
{
    for (size_t i = 0; i < manifest->n_fields; i++) {
        if (manifest->fields[i].id == id) {
            return manifest->fields[i].origin;
        }
    }

    return OUT_IF_PRESENT;
}

bool out_wanted(out_field_id id, const bool *wanted)
{
    return wanted[id];
}

/* Gives the values the widest row of the run occupies. A field left out is not measured,
 * so a run without the pairwise fields sizes its buffers from the per-base ones. */
size_t out_widest(size_t cap, const bool *wanted)
{
    size_t widest = 0;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        size_t width = wanted[id] ? out_values(id, cap, cap) : 0;

        widest = width > widest ? width : widest;
    }

    return widest;
}

/* Narrows one marker to the type a field is stored in. Refuses NaN in a type that has
 * none, where the cast is undefined. */
static int narrow_marker(out_field_id id, double marker, out_value *value)
{
    if (isnan(marker) && OUT_FIELDS[id].stored != OUT_F32) {
        return -1;
    }

    switch (OUT_FIELDS[id].stored) {
        case OUT_F32:      value->f32 = (float)marker;    return 0;
        case OUT_U64:      value->u64 = (uint64_t)marker; return 0;
        case OUT_I8:       value->i8  = (int8_t)marker;   return 0;
        case OUT_N_STORED: break;
    }

    return -1;
}

int out_fill_value(out_field_id id, out_value *value)
{
    return narrow_marker(id, OUT_FIELDS[id].fill, value);
}

bool out_values_needed(out_field_id id)
{
    return OUT_FIELDS[id].from_ref;
}

size_t out_stored_bytes(out_field_id id)
{
    switch (OUT_FIELDS[id].stored) {
        case OUT_F32:      return sizeof(float);
        case OUT_U64:      return sizeof(uint64_t);
        case OUT_I8:       return sizeof(int8_t);
        case OUT_N_STORED: break;
    }

    return 0;
}

size_t out_widest_bytes(void)
{
    size_t widest = 0;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        size_t bytes = out_stored_bytes(id);

        widest = bytes > widest ? bytes : widest;
    }

    return widest;
}

int out_dims(out_field_id id, int32_t n_refs, size_t cap, size_t *dims)
{
    shape_extents row     = OUT_FIELDS[id].row(cap, cap);
    int           extents = shape_rank(row);
    int           rank    = 0;

    if (OUT_FIELDS[id].per_ref) {
        dims[rank++] = (size_t)n_refs;
    }

    for (int i = 0; i < extents; i++) {
        dims[rank++] = row.dim[i];
    }

    return rank;
}

/* Gives the rank, discarding the dimensions an empty run produces. */
int out_rank(out_field_id id)
{
    size_t dims[OUT_RANK_MAX];

    return out_dims(id, 0, 0, dims);
}

/* ------------------------------------------------------------------------ */
/* The layout, described                                                     */
/* ------------------------------------------------------------------------ */

static const char *stored_name(out_stored stored)
{
    switch (stored) {
        case OUT_F32:      return "float32";
        case OUT_U64:      return "uint64";
        case OUT_I8:       return "int8";
        case OUT_N_STORED: break;
    }

    return "unknown";
}

/* Writes a field's fill into out, as the documentation names it. */
static void fill_name(out_field_id id, char *out, size_t len)
{
    double fill = OUT_FIELDS[id].fill;

    if (isnan(fill)) {
        snprintf(out, len, "nan");
    } else {
        snprintf(out, len, "%g", fill);
    }
}

/* Returns a sentence as JSON, or null where a field has none. Escapes the characters
 * that would otherwise end the string early. */
static void print_detail(FILE *out, const char *detail)
{
    if (!detail) {
        fputs("null", out);
        return;
    }

    fputc('"', out);

    for (const char *at = detail; *at; at++) {
        switch (*at) {
            case '"':  fputs("\\\"", out); break;
            case '\\': fputs("\\\\", out); break;
            case '\n': fputs("\\n", out);  break;
            default:   fputc(*at, out);    break;
        }
    }

    fputc('"', out);
}

/* Writes the dataset's extents as a JSON array of strings: "n" for the reference
 * dimension, then the row's extents, each as its symbol where it varies with the run and
 * as its number where it does not. */
static void print_extents(FILE *out, out_field_id id)
{
    shape_extents row  = OUT_FIELDS[id].row(0, 0);
    int           rank = shape_rank(row);

    fputc('[', out);

    if (OUT_FIELDS[id].per_ref) {
        fprintf(out, "\"n\"%s", rank > 0 ? ", " : "");
    }

    for (int i = 0; i < rank; i++) {
        const char *symbol = shape_symbol(OUT_FIELDS[id].row, i);

        if (symbol) {
            fprintf(out, "\"%s\"", symbol);
        } else {
            fprintf(out, "\"%zu\"", row.dim[i]);
        }

        if (i + 1 < rank) {
            fputs(", ", out);
        }
    }

    fputc(']', out);
}

/* Every attribute the output carries, as the objects of a JSON array. */
static void dump_attributes(FILE *out)
{
    for (out_attr_id id = 0; id < OUT_N_ATTRS; id++) {
        fprintf(out, "    {\n      \"name\": \"%s\",\n      \"detail\": ",
                OUT_ATTRIBUTES[id].name);

        print_detail(out, OUT_ATTRIBUTES[id].detail);

        fprintf(out, "\n    }%s\n", id + 1 < OUT_N_ATTRS ? "," : "");
    }
}

/* cmuts_version is the version of the program dumping this, which is not the version
 * attribute above: that one is written into a file and says what produced it. */
void out_dump_layout(FILE *out, const char *program, const out_manifest *manifest)
{
    fprintf(out, "{\n  \"program\": \"%s\",\n  \"cmuts_version\": \"%s\",\n"
                 "  \"attributes\": [\n", program, CMUTS_VERSION);

    dump_attributes(out);

    fprintf(out, "  ],\n  \"fields\": [\n");

    for (size_t i = 0; i < manifest->n_fields; i++) {
        out_field_id     id     = manifest->fields[i].id;
        const out_field *field  = &OUT_FIELDS[id];
        const char      *detail = manifest->fields[i].detail;
        const char      *needs  = manifest->fields[i].condition;
        char             fill[32];

        fill_name(id, fill, sizeof fill);

        fprintf(out,
                "    {\n"
                "      \"name\": \"%s\",\n"
                "      \"row\": \"%s\",\n"
                "      \"per_reference\": %s,\n"
                "      \"rank\": %d,\n"
                "      \"extents\": ",
                field->name, shape_name(field->row),
                field->per_ref ? "true" : "false", out_rank(id));

        print_extents(out, id);

        fprintf(out,
                ",\n"
                "      \"condition\": %s%s%s,\n"
                "      \"type\": \"%s\",\n"
                "      \"absent\": \"%s\",\n"
                "      \"detail\": ",
                needs ? "\"" : "", needs ? needs : "null", needs ? "\"" : "",
                stored_name(field->stored), fill);

        print_detail(out, detail);

        fprintf(out, "\n    }%s\n", i + 1 < manifest->n_fields ? "," : "");
    }

    fprintf(out, "  ]\n}\n");
}
