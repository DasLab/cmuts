/* cli.c -- command line parsing, driven by tables of options and positionals.
 *
 * Every argument the program accepts is described once, in OPTIONS or
 * POSITIONALS. The short option string, the getopt_long table, assignment,
 * bounds checking, which arguments are mandatory, the usage line, the help
 * text and the JSON description all derive from those tables, so nothing can
 * fall out of step with them and adding an argument means adding a row.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "cli.h"

#include <getopt.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* Long options without a short form need an identifier getopt_long can return;
 * anything above the character range will do. */
#define OPTION_ID_BASE 256

typedef enum {
    OPT_FLAG,    /* takes no argument; sets a bool */
    OPT_STRING,
    OPT_SIZE,
    OPT_INT,
} opt_type;

typedef struct {
    const char *group;     /* heading this option appears under in the help */
    const char *name;      /* long form */
    char        key;       /* short form, or 0 for none */
    opt_type    type;
    size_t      offset;    /* destination within cli_args */
    const char *metavar;   /* argument placeholder; NULL when it takes none */
    const char *help;      /* one line, for the help output */
    const char *detail;    /* paragraph for manual pages; NULL to reuse help */
    bool        required;
    long        minimum;   /* bounds for the numeric types */
    long        maximum;
    bool        hidden;    /* kept out of the help, still described by JSON */
} cli_option;

typedef struct {
    const char *name;      /* identifier for documentation and completions */
    const char *metavar;   /* placeholder in the usage line */
    const char *help;
    const char *detail;
    size_t      offset;
    bool        required;
} cli_positional;

static const cli_option OPTIONS[] = {
    { "Input and output", "fasta", 'f', OPT_STRING,
      offsetof(cli_args, pipeline.fasta_path), "FASTA",
      "reference sequences",
      "Reference sequences, in the same order as the alignment file's header. "
      "Read as a single forward pass and never indexed, so the file may be of "
      "any size.",
      true, 0, 0, false },
    { "Input and output", "output", 'o', OPT_STRING,
      offsetof(cli_args, pipeline.output_path), "HDF5",
      "write results to this file",
      "Results are written as one row per reference, with references that "
      "received no reads left as NaN.",
      true, 0, 0, false },

    { "Performance", "workers", 'j', OPT_SIZE,
      offsetof(cli_args, pipeline.workers), "N",
      "threads running the processing step",
      "Reads are taken from a shared pool, so a worker is free to cross "
      "reference boundaries and no thread idles waiting for a reference of "
      "its own.",
      false, 1, 1024, false },
    { "Performance", "decode-threads", 0, OPT_INT,
      offsetof(cli_args, pipeline.decode_threads), "N",
      "htslib threads for BGZF decompression",
      "Parallelizes inflation only; reading and record parsing stay "
      "sequential. Worth raising when the loader is the bottleneck, and "
      "pointless on small files.",
      false, 0, 64, false },
    { "Performance", "queue-capacity", 0, OPT_SIZE,
      offsetof(cli_args, pipeline.queue_capacity), "N",
      "reads in transit at once",
      "Bounds how far the loader may run ahead of the workers, and with it "
      "how much memory reads in flight occupy.",
      false, 1, 1 << 20, false },
    { "Performance", "batch", 0, OPT_SIZE,
      offsetof(cli_args, pipeline.batch), "N",
      "reads transferred per queue operation",
      "Larger batches amortize locking across more reads, at the cost of "
      "holding that many reads behind a slow one.",
      false, 1, 1 << 16, false },
    { "Performance", "live-refs", 0, OPT_SIZE,
      offsetof(cli_args, pipeline.live_refs), "N",
      "references in flight; 0 derives it from memory",
      "How far the loader may run ahead of a worker that stalls on one read. "
      "Zero derives a count from the longest reference and a memory budget.",
      false, 0, 1 << 16, false },

    { "Information", "help", 'h', OPT_FLAG,
      offsetof(cli_args, show_help), NULL,
      "show this help and exit", NULL, false, 0, 0, false },
    { "Information", "version", 'V', OPT_FLAG,
      offsetof(cli_args, show_version), NULL,
      "show the version and exit", NULL, false, 0, 0, false },
    { "Information", "dump-options", 0, OPT_FLAG,
      offsetof(cli_args, dump_options), NULL,
      "describe every argument as JSON and exit",
      "Intended for generating documentation and shell completions from the "
      "binary rather than from a separate description of it.",
      false, 0, 0, true },
};

static const cli_positional POSITIONALS[] = {
    { "alignment", "BAM", "coordinate-sorted alignments",
      "Read once, sequentially. Must be coordinate sorted, so that a "
      "reference is finished as soon as the reader moves past it.",
      offsetof(cli_args, pipeline.bam_path), true },
};

#define OPTION_COUNT     (sizeof OPTIONS / sizeof *OPTIONS)
#define POSITIONAL_COUNT (sizeof POSITIONALS / sizeof *POSITIONALS)

cli_args cli_defaults(void)
{
    return (cli_args){ .pipeline = pipeline_defaults() };
}

/* ------------------------------------------------------------------------ */
/* Table lookups                                                             */
/* ------------------------------------------------------------------------ */

static bool takes_argument(const cli_option *opt)
{
    return opt->type != OPT_FLAG;
}

static const cli_option *option_by_key(char key)
{
    for (size_t i = 0; i < OPTION_COUNT; i++)
        if (OPTIONS[i].key == key)
            return &OPTIONS[i];

    return NULL;
}

static const cli_option *option_by_id(int id)
{
    size_t index = (size_t)(id - OPTION_ID_BASE);

    return index < OPTION_COUNT ? &OPTIONS[index] : NULL;
}

/* ------------------------------------------------------------------------ */
/* getopt_long inputs, built from the table                                  */
/* ------------------------------------------------------------------------ */

static void build_short_options(char *out, size_t size)
{
    size_t n = 0;

    for (size_t i = 0; i < OPTION_COUNT && n + 3 < size; i++) {
        if (!OPTIONS[i].key)
            continue;

        out[n++] = OPTIONS[i].key;
        if (takes_argument(&OPTIONS[i]))
            out[n++] = ':';
    }

    out[n] = '\0';
}

static void build_long_options(struct option *out)
{
    for (size_t i = 0; i < OPTION_COUNT; i++) {
        out[i].name    = OPTIONS[i].name;
        out[i].has_arg = takes_argument(&OPTIONS[i]) ? required_argument : no_argument;
        out[i].flag    = NULL;
        out[i].val     = (int)(OPTION_ID_BASE + i);
    }

    out[OPTION_COUNT] = (struct option){ 0 };
}

/* ------------------------------------------------------------------------ */
/* Assignment                                                                */
/* ------------------------------------------------------------------------ */

static int parse_number(const cli_option *opt, const char *text, const char *program, long *out)
{
    char *end = NULL;
    long  n   = strtol(text, &end, 10);

    if (*text == '\0' || *end != '\0') {
        fprintf(stderr, "%s: --%s: \"%s\" is not a number\n", program, opt->name, text);
        return -1;
    }

    if (n < opt->minimum || n > opt->maximum) {
        fprintf(stderr, "%s: --%s: %ld is outside %ld..%ld\n",
                program, opt->name, n, opt->minimum, opt->maximum);
        return -1;
    }

    *out = n;
    return 0;
}

static int assign(const cli_option *opt, cli_args *args, const char *value, const char *program)
{
    char *field = (char *)args + opt->offset;
    long  number;

    switch (opt->type) {
        case OPT_FLAG:
            *(bool *)field = true;
            return 0;

        case OPT_STRING:
            *(const char **)field = value;
            return 0;

        case OPT_SIZE:
            if (parse_number(opt, value, program, &number) < 0)
                return -1;
            *(size_t *)field = (size_t)number;
            return 0;

        case OPT_INT:
            if (parse_number(opt, value, program, &number) < 0)
                return -1;
            *(int *)field = (int)number;
            return 0;
    }

    return -1;
}

/* ------------------------------------------------------------------------ */
/* Help                                                                      */
/* ------------------------------------------------------------------------ */

/* Renders an option's default from a defaults-initialised cli_args, so the
 * help can never advertise a value the program does not actually use. */
static void format_default(const cli_option *opt, const cli_args *defaults,
                           char *out, size_t size)
{
    const char *field = (const char *)defaults + opt->offset;

    switch (opt->type) {
        case OPT_SIZE:
            snprintf(out, size, " (default %zu)", *(const size_t *)field);
            break;
        case OPT_INT:
            snprintf(out, size, " (default %d)", *(const int *)field);
            break;
        default:
            out[0] = '\0';
            break;
    }
}

static void print_option(FILE *out, const cli_option *opt, const cli_args *defaults)
{
    char invocation[64];
    char suffix[48];
    int  n;

    if (opt->key)
        n = snprintf(invocation, sizeof invocation, "-%c, --%s", opt->key, opt->name);
    else
        n = snprintf(invocation, sizeof invocation, "    --%s", opt->name);

    if (opt->metavar && n > 0 && (size_t)n < sizeof invocation)
        snprintf(invocation + n, sizeof invocation - (size_t)n, " %s", opt->metavar);

    format_default(opt, defaults, suffix, sizeof suffix);
    fprintf(out, "  %-28s %s%s\n", invocation, opt->help, suffix);
}

/* Groups are collected by scanning rather than by assuming the table is sorted
 * by group, so a row added anywhere lands under the right heading. Hidden rows
 * do not count, so a group starts at its first visible option and a group with
 * none at all never prints a heading. */
static bool group_seen_before(size_t index)
{
    for (size_t i = 0; i < index; i++)
        if (!OPTIONS[i].hidden && strcmp(OPTIONS[i].group, OPTIONS[index].group) == 0)
            return true;

    return false;
}

static void print_usage_line(FILE *out, const char *program)
{
    fprintf(out, "usage: %s [options]", program);

    for (size_t i = 0; i < OPTION_COUNT; i++) {
        if (!OPTIONS[i].required)
            continue;

        if (OPTIONS[i].key)
            fprintf(out, " -%c %s", OPTIONS[i].key, OPTIONS[i].metavar);
        else
            fprintf(out, " --%s %s", OPTIONS[i].name, OPTIONS[i].metavar);
    }

    for (size_t i = 0; i < POSITIONAL_COUNT; i++)
        fprintf(out, POSITIONALS[i].required ? " %s" : " [%s]", POSITIONALS[i].metavar);

    fputc('\n', out);
}

static void print_positionals(FILE *out)
{
    fprintf(out, "\nArguments\n");

    for (size_t i = 0; i < POSITIONAL_COUNT; i++)
        fprintf(out, "  %-28s %s\n", POSITIONALS[i].metavar, POSITIONALS[i].help);
}

void cli_usage(FILE *out, const char *program)
{
    cli_args defaults = cli_defaults();

    fprintf(out, "cmuts %s -- per-reference accumulation over a BAM file\n\n", CMUTS_VERSION);
    print_usage_line(out, program);
    print_positionals(out);

    for (size_t i = 0; i < OPTION_COUNT; i++) {
        if (OPTIONS[i].hidden || group_seen_before(i))
            continue;

        fprintf(out, "\n%s\n", OPTIONS[i].group);

        for (size_t j = i; j < OPTION_COUNT; j++)
            if (!OPTIONS[j].hidden && strcmp(OPTIONS[j].group, OPTIONS[i].group) == 0)
                print_option(out, &OPTIONS[j], &defaults);
    }
}

/* ------------------------------------------------------------------------ */
/* JSON description                                                          */
/* ------------------------------------------------------------------------ */

static void print_json_string(FILE *out, const char *text)
{
    if (!text) {
        fputs("null", out);
        return;
    }

    fputc('"', out);
    for (const unsigned char *p = (const unsigned char *)text; *p; p++) {
        switch (*p) {
            case '"':  fputs("\\\"", out); break;
            case '\\': fputs("\\\\", out); break;
            case '\n': fputs("\\n", out);  break;
            case '\t': fputs("\\t", out);  break;
            default:
                if (*p < 0x20)
                    fprintf(out, "\\u%04x", *p);
                else
                    fputc(*p, out);
        }
    }
    fputc('"', out);
}

static const char *type_name(opt_type type)
{
    switch (type) {
        case OPT_FLAG:   return "flag";
        case OPT_STRING: return "string";
        case OPT_SIZE:   return "size";
        case OPT_INT:    return "int";
    }

    return "unknown";
}

static void print_json_default(FILE *out, const cli_option *opt, const cli_args *defaults)
{
    const char *field = (const char *)defaults + opt->offset;

    switch (opt->type) {
        case OPT_FLAG:   fputs(*(const bool *)field ? "true" : "false", out); break;
        case OPT_STRING: print_json_string(out, *(const char *const *)field); break;
        case OPT_SIZE:   fprintf(out, "%zu", *(const size_t *)field);         break;
        case OPT_INT:    fprintf(out, "%d", *(const int *)field);             break;
    }
}

static void print_json_bounds(FILE *out, const cli_option *opt)
{
    if (opt->type != OPT_SIZE && opt->type != OPT_INT) {
        fputs("      \"minimum\": null, \"maximum\": null\n", out);
        return;
    }

    fprintf(out, "      \"minimum\": %ld, \"maximum\": %ld\n", opt->minimum, opt->maximum);
}

static void print_json_option(FILE *out, const cli_option *opt, const cli_args *defaults,
                              bool last)
{
    char key[2] = { opt->key, '\0' };

    fputs("    {\n", out);
    fputs("      \"name\": ", out);         print_json_string(out, opt->name);
    fputs(",\n      \"short\": ", out);     print_json_string(out, opt->key ? key : NULL);
    fputs(",\n      \"group\": ", out);     print_json_string(out, opt->group);
    fprintf(out, ",\n      \"type\": \"%s\"", type_name(opt->type));
    fputs(",\n      \"metavar\": ", out);   print_json_string(out, opt->metavar);
    fputs(",\n      \"help\": ", out);      print_json_string(out, opt->help);
    fputs(",\n      \"detail\": ", out);    print_json_string(out, opt->detail);
    fprintf(out, ",\n      \"required\": %s", opt->required ? "true" : "false");
    fprintf(out, ",\n      \"hidden\": %s", opt->hidden ? "true" : "false");
    fputs(",\n      \"default\": ", out);   print_json_default(out, opt, defaults);
    fputs(",\n", out);                      print_json_bounds(out, opt);
    fprintf(out, "    }%s\n", last ? "" : ",");
}

static void print_json_positional(FILE *out, const cli_positional *pos, bool last)
{
    fputs("    {\n", out);
    fputs("      \"name\": ", out);       print_json_string(out, pos->name);
    fputs(",\n      \"metavar\": ", out); print_json_string(out, pos->metavar);
    fputs(",\n      \"help\": ", out);    print_json_string(out, pos->help);
    fputs(",\n      \"detail\": ", out);  print_json_string(out, pos->detail);
    fprintf(out, ",\n      \"required\": %s\n", pos->required ? "true" : "false");
    fprintf(out, "    }%s\n", last ? "" : ",");
}

void cli_dump_options(FILE *out)
{
    cli_args defaults = cli_defaults();

    fputs("{\n", out);
    fputs("  \"program\": \"cmuts\",\n", out);
    fprintf(out, "  \"version\": \"%s\",\n", CMUTS_VERSION);

    fputs("  \"options\": [\n", out);
    for (size_t i = 0; i < OPTION_COUNT; i++)
        print_json_option(out, &OPTIONS[i], &defaults, i + 1 == OPTION_COUNT);
    fputs("  ],\n", out);

    fputs("  \"positionals\": [\n", out);
    for (size_t i = 0; i < POSITIONAL_COUNT; i++)
        print_json_positional(out, &POSITIONALS[i], i + 1 == POSITIONAL_COUNT);
    fputs("  ]\n", out);

    fputs("}\n", out);
}

/* ------------------------------------------------------------------------ */
/* Parsing                                                                   */
/* ------------------------------------------------------------------------ */

static cli_status check_required_options(const bool *seen, const char *program)
{
    for (size_t i = 0; i < OPTION_COUNT; i++) {
        if (!OPTIONS[i].required || seen[i])
            continue;

        fprintf(stderr, "%s: missing required option --%s (%s)\n",
                program, OPTIONS[i].name, OPTIONS[i].help);
        return CLI_ERROR;
    }

    return CLI_OK;
}

static void print_expected_positionals(FILE *out)
{
    for (size_t i = 0; i < POSITIONAL_COUNT; i++)
        fprintf(out, "%s%s", i ? " " : "", POSITIONALS[i].metavar);
}

static cli_status take_positionals(int argc, char **argv, cli_args *args, const char *program)
{
    int given = argc - optind;

    if (given != (int)POSITIONAL_COUNT) {
        fprintf(stderr, "%s: expected ", program);
        print_expected_positionals(stderr);
        fprintf(stderr, ", got %d argument%s\n", given, given == 1 ? "" : "s");
        return CLI_ERROR;
    }

    for (size_t i = 0; i < POSITIONAL_COUNT; i++) {
        const char **field = (const char **)((char *)args + POSITIONALS[i].offset);

        *field = argv[optind + (int)i];
    }

    return CLI_OK;
}

/* Requests that are answered entirely by the tables, before any argument is
 * required to be present. */
static bool answered_immediately(const cli_args *args, const char *program)
{
    if (args->show_help) {
        cli_usage(stdout, program);
        return true;
    }

    if (args->show_version) {
        printf("cmuts %s\n", CMUTS_VERSION);
        return true;
    }

    if (args->dump_options) {
        cli_dump_options(stdout);
        return true;
    }

    return false;
}

cli_status cli_parse(int argc, char **argv, cli_args *args)
{
    struct option longopts[OPTION_COUNT + 1];
    char          shortopts[2 * OPTION_COUNT + 1];
    bool          seen[OPTION_COUNT];
    const char   *program = argv[0];
    int           found;

    *args = cli_defaults();
    memset(seen, 0, sizeof seen);
    build_short_options(shortopts, sizeof shortopts);
    build_long_options(longopts);

    opterr = 0;
    optind = 1;

    while ((found = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1) {
        const cli_option *opt = found < OPTION_ID_BASE ? option_by_key((char)found)
                                                       : option_by_id(found);

        if (!opt) {
            fprintf(stderr, "%s: unrecognized option; try --help\n", program);
            return CLI_ERROR;
        }

        if (assign(opt, args, optarg, program) < 0)
            return CLI_ERROR;

        seen[opt - OPTIONS] = true;
    }

    if (answered_immediately(args, program))
        return CLI_DONE;

    if (check_required_options(seen, program) != CLI_OK)
        return CLI_ERROR;

    return take_positionals(argc, argv, args, program);
}
