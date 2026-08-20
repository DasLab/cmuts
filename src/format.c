/* output.c -- the output field table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "format.h"

#include <math.h>

#include "version.h"

const fmt_field FMT_FIELDS[FMT_N_FIELDS] = {
    [FMT_COVERAGE] = {
        .name    = "coverage",
        .detail  = "The number of reads in which this base was present.",
        .row     = shape_per_base,
        .per_ref = true,
        .stored  = FMT_F32,
        .fill    = 0.0,
    },
    [FMT_REACTIVITY] = {
        .name    = "reactivity",
        .detail  = "The reactivity at each base.",
        .row     = shape_per_base,
        .per_ref = true,
        .stored  = FMT_F32,
        .fill    = (double)NAN,
    },
    [FMT_ERROR] = {
        .name    = "error",
        .detail  = "Standard error of the reactivity values.",
        .row     = shape_per_base,
        .per_ref = true,
        .stored  = FMT_F32,
        .fill    = (double)NAN,
    },
    [FMT_LENGTHS] = {
        .name    = "reads/lengths",
        .detail  = "The number of reads passing all filters, binned by length.",
        .row     = shape_per_length,
        .per_ref = true,
        .stored  = FMT_U64,
        .fill    = 0.0,
    },
    [FMT_READS] = {
        .name    = "reads/counted",
        .detail  = "The number of reads passing all filters.",
        .row     = shape_none,
        .per_ref = true,
        .stored  = FMT_U64,
        .fill    = 0.0,
    },
    [FMT_REJECTED] = {
        .name    = "reads/rejected",
        .detail  = "The number of reads rejected by at least one filter, or which couldn't be modelled by the HMM.",
        .row     = shape_none,
        .per_ref = true,
        .stored  = FMT_U64,
        .fill    = 0.0,
    },
    [FMT_PAIRWISE_CORRELATION] = {
        .name     = "pairwise/correlation",
        .detail   = "The Pearson correlation of mutations between this pair of bases.",
        .row      = shape_per_pair,
        .per_ref  = true,
        .stored   = FMT_F32,
        .fill     = (double)NAN,
    },
    [FMT_PAIRWISE_CONDITIONAL] = {
        .name     = "pairwise/conditional",
        .detail   = "The probability that the base on the first axis was mutated in a read, given that the base on the second axis was.",
        .row      = shape_per_pair,
        .per_ref  = true,
        .stored   = FMT_F32,
        .fill     = (double)NAN,
    },
    [FMT_PAIRWISE_COVERAGE] = {
        .name     = "pairwise/coverage",
        .detail   = "The number of reads in which this pair of bases was present.",
        .row      = shape_per_pair,
        .per_ref  = true,
        .stored   = FMT_F32,
        .fill     = 0.0,
    },
    [FMT_NORM] = {
        .name     = "norm",
        .detail   = "The norm every rate in this file was divided by.",
        .row      = shape_none,
        .per_ref  = false,
        .stored   = FMT_F32,
        .fill     = (double)NAN,
    },
    [FMT_SEQUENCE] = {
        .name     = "sequence",
        .detail   = "The reference sequence: 0 for A, 1 for C, 2 for G, 3 for T, and -1 for any other base and for every column past the reference's end.",
        .row      = shape_per_base,
        .per_ref  = true,
        .from_ref = true,
        .stored   = FMT_I8,
        .fill     = -1.0,
    },
    [FMT_UNMAPPED] = {
        .name    = "reads/unmapped",
        .detail  = "The number of reads not aligned to any reference.",
        .row     = shape_none,
        .per_ref = false,
        .stored  = FMT_U64,
        .fill    = 0.0,
    },
};

const fmt_attribute FMT_ATTRIBUTES[FMT_N_ATTRS] = {
    [FMT_ATTR_PROGRAM] = {
        .name   = "program",
        .detail = "The name of the program that produced this file.",
    },
    [FMT_ATTR_VERSION] = {
        .name   = "version",
        .detail = "The version of `cmuts` that produced this file.",
    },
};

size_t fmt_values(fmt_field_id id, size_t len, size_t cap)
{
    return shape_values(FMT_FIELDS[id].row, len, cap);
}

void fmt_selection(const fmt_manifest *manifest, bool *wanted)
{
    for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
        wanted[id] = false;
    }

    for (size_t i = 0; i < manifest->n_fields; i++) {
        wanted[manifest->fields[i].id] = true;
    }
}

/* Adds one field to the reads, keeping each field a single entry. */
static void reads_add(fmt_reads *reads, fmt_field_id id, bool required)
{
    for (size_t i = 0; i < reads->n_fields; i++) {
        if (reads->fields[i].id == id) {
            if (required) {
                reads->fields[i].required = true;
            }
            return;
        }
    }

    reads->fields[reads->n_fields++] = (fmt_read){ id, required };
}

void fmt_reads_of(const fmt_manifest *manifest, fmt_reads *reads)
{
    reads->n_fields = 0;

    for (size_t i = 0; i < manifest->n_fields; i++) {
        const fmt_written *field = &manifest->fields[i];

        for (const fmt_field_id *dep = field->depends;
             dep && *dep != FMT_N_FIELDS; dep++) {
            reads_add(reads, *dep, field->required);
        }
    }
}

void fmt_reads_selection(const fmt_reads *reads, bool *wanted)
{
    for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
        wanted[id] = false;
    }

    for (size_t i = 0; i < reads->n_fields; i++) {
        wanted[reads->fields[i].id] = true;
    }
}

bool fmt_read_required(const fmt_reads *reads, fmt_field_id id)
{
    for (size_t i = 0; i < reads->n_fields; i++) {
        if (reads->fields[i].id == id) {
            return reads->fields[i].required;
        }
    }

    return false;
}

bool fmt_wanted(fmt_field_id id, const bool *wanted)
{
    return wanted[id];
}

/* Gives the values the widest row of the run occupies. A field left out is not measured,
 * so a run without the pairwise fields sizes its buffers from the per-base ones. */
size_t fmt_widest(size_t cap, const bool *wanted)
{
    size_t widest = 0;

    for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
        size_t width = wanted[id] ? fmt_values(id, cap, cap) : 0;

        widest = width > widest ? width : widest;
    }

    return widest;
}

/* Narrows one marker to the type a field is stored in. Refuses NaN in a type that has
 * none, where the cast is undefined. */
static int narrow_marker(fmt_field_id id, double marker, fmt_value *value)
{
    if (isnan(marker) && FMT_FIELDS[id].stored != FMT_F32) {
        return -1;
    }

    switch (FMT_FIELDS[id].stored) {
        case FMT_F32:      value->f32 = (float)marker;    return 0;
        case FMT_U64:      value->u64 = (uint64_t)marker; return 0;
        case FMT_I8:       value->i8  = (int8_t)marker;   return 0;
        case FMT_N_STORED: break;
    }

    return -1;
}

int fmt_fill_value(fmt_field_id id, fmt_value *value)
{
    return narrow_marker(id, FMT_FIELDS[id].fill, value);
}

bool fmt_values_needed(fmt_field_id id)
{
    return FMT_FIELDS[id].from_ref;
}

size_t fmt_stored_bytes(fmt_field_id id)
{
    switch (FMT_FIELDS[id].stored) {
        case FMT_F32:      return sizeof(float);
        case FMT_U64:      return sizeof(uint64_t);
        case FMT_I8:       return sizeof(int8_t);
        case FMT_N_STORED: break;
    }

    return 0;
}

size_t fmt_widest_bytes(void)
{
    size_t widest = 0;

    for (fmt_field_id id = 0; id < FMT_N_FIELDS; id++) {
        size_t bytes = fmt_stored_bytes(id);

        widest = bytes > widest ? bytes : widest;
    }

    return widest;
}

int fmt_dims(fmt_field_id id, int32_t n_refs, size_t cap, size_t *dims)
{
    shape_extents row     = FMT_FIELDS[id].row(cap, cap);
    int           extents = shape_rank(row);
    int           rank    = 0;

    if (FMT_FIELDS[id].per_ref) {
        dims[rank++] = (size_t)n_refs;
    }

    for (int i = 0; i < extents; i++) {
        dims[rank++] = row.dim[i];
    }

    return rank;
}

/* Gives the rank, discarding the dimensions an empty run produces. */
int fmt_rank(fmt_field_id id)
{
    size_t dims[FMT_RANK_MAX];

    return fmt_dims(id, 0, 0, dims);
}

/* ------------------------------------------------------------------------ */
/* The layout, described                                                     */
/* ------------------------------------------------------------------------ */

static const char *stored_name(fmt_stored stored)
{
    switch (stored) {
        case FMT_F32:      return "float32";
        case FMT_U64:      return "uint64";
        case FMT_I8:       return "int8";
        case FMT_N_STORED: break;
    }

    return "unknown";
}

/* Writes a field's fill into out, as the documentation names it. */
static void fill_name(fmt_field_id id, char *out, size_t len)
{
    double fill = FMT_FIELDS[id].fill;

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
static void print_extents(FILE *out, fmt_field_id id)
{
    shape_extents row  = FMT_FIELDS[id].row(0, 0);
    int           rank = shape_rank(row);

    fputc('[', out);

    if (FMT_FIELDS[id].per_ref) {
        fprintf(out, "\"n\"%s", rank > 0 ? ", " : "");
    }

    for (int i = 0; i < rank; i++) {
        const char *symbol = shape_symbol(FMT_FIELDS[id].row, i);

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
    for (fmt_attr_id id = 0; id < FMT_N_ATTRS; id++) {
        fprintf(out, "    {\n      \"name\": \"%s\",\n      \"detail\": ",
                FMT_ATTRIBUTES[id].name);

        print_detail(out, FMT_ATTRIBUTES[id].detail);

        fprintf(out, "\n    }%s\n", id + 1 < FMT_N_ATTRS ? "," : "");
    }
}

/* cmuts_version is the version of the program dumping this, which is not the version
 * attribute above: that one is written into a file and says what produced it. */
void fmt_dump_layout(FILE *out, const char *program, const fmt_manifest *manifest)
{
    fprintf(out, "{\n  \"program\": \"%s\",\n  \"cmuts_version\": \"%s\",\n"
                 "  \"attributes\": [\n", program, CMUTS_VERSION);

    dump_attributes(out);

    fprintf(out, "  ],\n  \"fields\": [\n");

    for (size_t i = 0; i < manifest->n_fields; i++) {
        fmt_field_id     id    = manifest->fields[i].id;
        const fmt_field *field = &FMT_FIELDS[id];
        const char      *note  = manifest->fields[i].note;
        const char      *needs = manifest->fields[i].condition;
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
                field->per_ref ? "true" : "false", fmt_rank(id));

        print_extents(out, id);

        fprintf(out,
                ",\n"
                "      \"condition\": %s%s%s,\n"
                "      \"type\": \"%s\",\n"
                "      \"absent\": \"%s\",\n"
                "      \"detail\": ",
                needs ? "\"" : "", needs ? needs : "null", needs ? "\"" : "",
                stored_name(field->stored), fill);

        print_detail(out, field->detail);

        fprintf(out, ",\n      \"note\": ");

        print_detail(out, note);

        fprintf(out, "\n    }%s\n", i + 1 < manifest->n_fields ? "," : "");
    }

    fprintf(out, "  ]\n}\n");
}
