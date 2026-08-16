/* output.c -- the output field table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "output.h"

#include "version.h"

const out_field OUT_FIELDS[OUT_N_FIELDS] = {
    [OUT_COVERAGE] = {
        .name    = "coverage",
        .row     = shape_per_base,
        .per_ref = true,
        .stored  = OUT_F32,
        .absent  = OUT_ZERO,
        .detail  = "The number of reads in which this base was present, weighted by PHRED scores.",
    },
    [OUT_REACTIVITY] = {
        .name    = "reactivity",
        .row     = shape_per_base,
        .per_ref = true,
        .stored  = OUT_F32,
        .absent  = OUT_NAN,
        .detail  = "The mutation rate at this base, weighted by PHRED scores and in accordance with the HMM parameters.",
    },
    [OUT_ERROR] = {
        .name    = "error",
        .row     = shape_per_base,
        .per_ref = true,
        .stored  = OUT_F32,
        .absent  = OUT_NAN,
        .detail  = "Standard error of the reactivity values. Purely the statistical error introduced by finite read depths; does not account for experimental or systemic errors.",
    },
    [OUT_LENGTHS] = {
        .name    = "reads/lengths",
        .row     = shape_per_length,
        .per_ref = true,
        .stored  = OUT_U64,
        .absent  = OUT_ZERO,
        .detail  = "The number of reads passing all filters, binned by length.",
    },
    [OUT_READS] = {
        .name    = "reads/counted",
        .row     = shape_none,
        .per_ref = true,
        .stored  = OUT_U64,
        .absent  = OUT_ZERO,
        .detail  = "The number of reads passing all filters.",
    },
    [OUT_REJECTED] = {
        .name    = "reads/rejected",
        .row     = shape_none,
        .per_ref = true,
        .stored  = OUT_U64,
        .absent  = OUT_ZERO,
        .detail  = "The number of reads that contributed nothing: rejected by at least one filter, or carrying something the pair HMM's rates give no alignment of.",
    },
    [OUT_PAIRWISE_CORRELATION] = {
        .name     = "pairwise/correlation",
        .row      = shape_per_pair,
        .per_ref  = true,
        .stored   = OUT_F32,
        .absent   = OUT_NAN,
        .optional = true,
        .detail   = "The correlation between two positions being modified in the same read, as the Pearson coefficient of the two binary variables. NaN where the reads are too few, and where either position is modified in all of them or in none. The diagonal is a position against itself, which falls short of one by however much of its variance the base calls leave unsettled; divide a correlation by the square root of the two diagonals to take that out.",
    },
    [OUT_PAIRWISE_COVERAGE] = {
        .name     = "pairwise/coverage",
        .row      = shape_per_pair,
        .per_ref  = true,
        .stored   = OUT_F32,
        .absent   = OUT_ZERO,
        .optional = true,
        .detail   = "The evidence behind each correlation: the reads reaching both positions.",
    },
    [OUT_UNMAPPED] = {
        .name    = "reads/unmapped",
        .row     = shape_none,
        .per_ref = false,
        .stored  = OUT_U64,
        .absent  = OUT_ZERO,
        .detail  = "The number of reads not aligned to any reference.",
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

const bool OUT_REQUIRED_ONLY[OUT_N_FIELDS] = { false };

bool out_wanted(out_field_id id, const bool *wanted)
{
    return !OUT_FIELDS[id].optional || !wanted || wanted[id];
}

/* Gives the values the widest row of the run occupies. A field left out is not measured,
 * so a run without the pairwise fields sizes its buffers from the per-base ones and pays
 * nothing for the square it did not ask for. */
size_t out_widest(size_t cap, const bool *wanted)
{
    size_t widest = 0;

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        size_t width = out_wanted(id, wanted) ? out_values(id, cap, cap) : 0;

        widest = width > widest ? width : widest;
    }

    return widest;
}

size_t out_stored_bytes(out_field_id id)
{
    switch (OUT_FIELDS[id].stored) {
        case OUT_F32:      return sizeof(float);
        case OUT_U64:      return sizeof(uint64_t);
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
        case OUT_N_STORED: break;
    }

    return "unknown";
}

/* Names the value read back where the run wrote nothing, which is not the same as the
 * type it is stored in: zero for a count, NaN for a rate. */
static const char *absent_name(out_absent absent)
{
    switch (absent) {
        case OUT_ZERO:     return "zero";
        case OUT_NAN:      return "nan";
        case OUT_N_ABSENT: break;
    }

    return "unknown";
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
void out_dump_layout(FILE *out)
{
    fprintf(out, "{\n  \"cmuts_version\": \"%s\",\n  \"attributes\": [\n", CMUTS_VERSION);

    dump_attributes(out);

    fprintf(out, "  ],\n  \"fields\": [\n");

    for (out_field_id id = 0; id < OUT_N_FIELDS; id++) {
        const out_field *field = &OUT_FIELDS[id];

        fprintf(out,
                "    {\n"
                "      \"name\": \"%s\",\n"
                "      \"row\": \"%s\",\n"
                "      \"per_reference\": %s,\n"
                "      \"rank\": %d,\n"
                "      \"type\": \"%s\",\n"
                "      \"absent\": \"%s\",\n"
                "      \"detail\": ",
                field->name, shape_name(field->row),
                field->per_ref ? "true" : "false", out_rank(id),
                stored_name(field->stored), absent_name(field->absent));

        print_detail(out, field->detail);

        fprintf(out, "\n    }%s\n", id + 1 < OUT_N_FIELDS ? "," : "");
    }

    fprintf(out, "  ]\n}\n");
}
