/* cli.c -- command line parsing, driven by a program's option table.
 *
 * Every argument a program accepts is described once, in its cli_spec. The
 * short option string, the getopt_long table, assignment, bounds checking,
 * which arguments are mandatory, the usage line, the help text and the JSON
 * description all derive from it, so nothing can fall out of step with it and
 * adding an argument means adding a row.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "cli.h"

#include <float.h>
#include <getopt.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Long options without a short form need an identifier getopt_long can return;
 * anything above the character range will do. */
#define OPTION_ID_BASE 256

/* What getopt returns for an option given without the value it needs. */
#define MISSING_VALUE ':'

/* Enough for a long option as anyone would write one. */
#define REJECTED_MAX 64

/* Enough for a placeholder and the ellipsis a variadic one carries. */
#define METAVAR_MAX 32

/* Enough for both forms of an option and the placeholder that follows them. */
#define INVOCATION_MAX 64

/* Enough for the note naming an option's default. A string default is written
 * out in full, so this is the one of these a row could outgrow, and the note is
 * cut short where it does rather than anything worse. */
#define DEFAULT_NOTE_MAX 64

/* Enough for the largest value an option will take, written out. */
#define CEILING_MAX 32

/* The column an option's help begins in. The arguments are listed against the
 * same one, and a difference between the two would read as a step partway down
 * the help rather than as anything deliberate. */
#define HELP_COLUMN 28

/* getopt records a short option in optopt, and for a long one leaves the
 * identifier it was given instead, which names nothing a reader would know. A
 * long one is read back from the word getopt stopped on. */
static void report_rejected_option(const cli_spec *spec, int argc, char **argv,
                                   bool needs_value)
{
    char rejected[REJECTED_MAX];

    if (optopt > 0 && optopt < OPTION_ID_BASE)
        snprintf(rejected, sizeof rejected, "-%c", optopt);
    else if (optind >= 1 && optind <= argc)
        snprintf(rejected, sizeof rejected, "%s", argv[optind - 1]);
    else
        snprintf(rejected, sizeof rejected, "an option");

    if (needs_value)
        fprintf(stderr, "%s: %s needs a value\n", spec->program, rejected);
    else
        fprintf(stderr, "%s: unrecognized option %s\n", spec->program, rejected);
}

static bool takes_argument(const cli_option *opt)
{
    return opt->type != OPT_FLAG;
}

/* An option that answers and exits has no destination, so nothing is written
 * for it and its offset is never used. */
static bool stores_a_value(const cli_option *opt)
{
    return opt->action == CLI_STORE;
}

/* ------------------------------------------------------------------------ */
/* Table lookups                                                             */
/* ------------------------------------------------------------------------ */

static const cli_option *option_by_key(const cli_spec *spec, char key)
{
    for (size_t i = 0; i < spec->n_options; i++)
        if (spec->options[i].key == key)
            return &spec->options[i];

    return NULL;
}

static const cli_option *option_by_id(const cli_spec *spec, int id)
{
    size_t index = (size_t)(id - OPTION_ID_BASE);

    return index < spec->n_options ? &spec->options[index] : NULL;
}

/* ------------------------------------------------------------------------ */
/* getopt_long inputs, built from the table                                  */
/* ------------------------------------------------------------------------ */

/* The leading colon has getopt report a missing value apart from an option it
 * does not know, which are different mistakes and want different words. */
static void build_short_options(const cli_spec *spec, char *out, size_t size)
{
    size_t n = 0;

    out[n++] = ':';

    for (size_t i = 0; i < spec->n_options && n + 3 < size; i++) {
        if (!spec->options[i].key)
            continue;

        out[n++] = spec->options[i].key;
        if (takes_argument(&spec->options[i]))
            out[n++] = ':';
    }

    out[n] = '\0';
}

static void build_long_options(const cli_spec *spec, struct option *out)
{
    for (size_t i = 0; i < spec->n_options; i++) {
        out[i].name    = spec->options[i].name;
        out[i].has_arg = takes_argument(&spec->options[i]) ? required_argument
                                                           : no_argument;
        out[i].flag    = NULL;
        out[i].val     = (int)(OPTION_ID_BASE + i);
    }

    out[spec->n_options] = (struct option){ 0 };
}

/* ------------------------------------------------------------------------ */
/* Assignment                                                                */
/* ------------------------------------------------------------------------ */

static void print_choices(FILE *out, const cli_option *opt)
{
    for (const cli_choice *choice = opt->choices; choice->name; choice++)
        fprintf(out, "%s%s", choice == opt->choices ? "" : "|", choice->name);
}

static int parse_choice(const cli_option *opt, const char *text, const char *program,
                        int *out)
{
    for (const cli_choice *choice = opt->choices; choice->name; choice++)
        if (strcmp(choice->name, text) == 0) {
            *out = choice->value;
            return 0;
        }

    fprintf(stderr, "%s: --%s: \"%s\" is not one of ", program, opt->name, text);
    print_choices(stderr, opt);
    fputc('\n', stderr);

    return -1;
}

/* The largest whole number an option will take. One with no ceiling of its own
 * is still held to what its destination can carry, so that no value is accepted
 * that the field would silently truncate. */
static long integer_ceiling(const cli_option *opt)
{
    if (opt->maximum != CLI_UNBOUNDED)
        return opt->maximum;

    return opt->type == OPT_INT ? INT32_MAX : LONG_MAX;
}

static double real_ceiling(const cli_option *opt)
{
    return opt->maximum != CLI_UNBOUNDED ? (double)opt->maximum : DBL_MAX;
}

/* The largest value an option will take, written out. */
static const char *ceiling_text(const cli_option *opt, char *out, size_t size)
{
    if (opt->type == OPT_DOUBLE)
        snprintf(out, size, "%g", real_ceiling(opt));
    else
        snprintf(out, size, "%ld", integer_ceiling(opt));

    return out;
}

/* A range is quoted only where the row set one. An option with no ceiling of
 * its own has a floor to fall short of and a destination to outgrow, and each
 * is one bound rather than a range: a value refused for being negative has no
 * business hearing what the widest one is. The text the caller wrote is echoed
 * rather than the number read out of it, so a value comes back as given. */
static void report_out_of_range(const cli_option *opt, const char *text,
                                const char *program, bool below)
{
    char limit[CEILING_MAX];

    if (opt->maximum != CLI_UNBOUNDED)
        fprintf(stderr, "%s: --%s: %s is outside %ld..%ld\n",
                program, opt->name, text, opt->minimum, opt->maximum);
    else if (below)
        fprintf(stderr, "%s: --%s: %s is below %ld\n",
                program, opt->name, text, opt->minimum);
    else
        fprintf(stderr, "%s: --%s: %s is above %s\n",
                program, opt->name, text, ceiling_text(opt, limit, sizeof limit));
}

static int parse_number(const cli_option *opt, const char *text, const char *program,
                        long *out)
{
    char *end = NULL;
    long  n   = strtol(text, &end, 10);

    if (*text == '\0' || *end != '\0') {
        fprintf(stderr, "%s: --%s: \"%s\" is not a number\n", program, opt->name, text);
        return -1;
    }

    if (n < opt->minimum || n > integer_ceiling(opt)) {
        report_out_of_range(opt, text, program, n < opt->minimum);
        return -1;
    }

    *out = n;
    return 0;
}

static int parse_double(const cli_option *opt, const char *text, const char *program,
                        double *out)
{
    char  *end = NULL;
    double n   = strtod(text, &end);

    if (*text == '\0' || *end != '\0') {
        fprintf(stderr, "%s: --%s: \"%s\" is not a number\n", program, opt->name, text);
        return -1;
    }

    if (n < (double)opt->minimum || n > real_ceiling(opt)) {
        report_out_of_range(opt, text, program, n < (double)opt->minimum);
        return -1;
    }

    *out = n;
    return 0;
}

static int assign(const cli_option *opt, void *args, const char *value,
                  const char *program)
{
    char  *field = (char *)args + opt->offset;
    long   number;
    double real;
    int    choice;

    switch (opt->type) {
        case OPT_DOUBLE:
            if (parse_double(opt, value, program, &real) < 0)
                return -1;
            *(double *)field = real;
            return 0;

        case OPT_ENUM:
            if (parse_choice(opt, value, program, &choice) < 0)
                return -1;
            *(int *)field = choice;
            return 0;

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

/* Renders an option's default from the spec's defaults, so the help can never
 * advertise a value the program does not actually use. */
static void format_default(const cli_option *opt, const void *defaults,
                           char *out, size_t size)
{
    const char *field = (const char *)defaults + opt->offset;

    out[0] = '\0';

    if (!stores_a_value(opt))
        return;

    /* An option that need not be applied says so rather than showing the value
     * that stands for not applying it. */
    if (opt->unset_label) {
        snprintf(out, size, " (default: %s)", opt->unset_label);
        return;
    }

    switch (opt->type) {
        case OPT_SIZE:
            snprintf(out, size, " (default %zu)", *(const size_t *)field);
            break;
        case OPT_INT:
            snprintf(out, size, " (default %d)", *(const int *)field);
            break;
        case OPT_DOUBLE:
            snprintf(out, size, " (default %g)", *(const double *)field);
            break;
        case OPT_STRING:
            if (*(const char *const *)field)
                snprintf(out, size, " (default %s)", *(const char *const *)field);
            break;
        case OPT_ENUM:
            for (const cli_choice *choice = opt->choices; choice->name; choice++)
                if (choice->value == *(const int *)field)
                    snprintf(out, size, " (default %s)", choice->name);
            break;
        default:
            break;
    }
}

static void print_option(FILE *out, const cli_option *opt, const void *defaults)
{
    char invocation[INVOCATION_MAX];
    char suffix[DEFAULT_NOTE_MAX];
    int  n;

    if (opt->key)
        n = snprintf(invocation, sizeof invocation, "-%c, --%s", opt->key, opt->name);
    else
        n = snprintf(invocation, sizeof invocation, "    --%s", opt->name);

    if (opt->metavar && n > 0 && (size_t)n < sizeof invocation)
        snprintf(invocation + n, sizeof invocation - (size_t)n, " %s", opt->metavar);

    fprintf(out, "  %-*s %s", HELP_COLUMN, invocation, opt->help);

    if (opt->choices) {
        fputs(" (", out);
        print_choices(out, opt);
        fputc(')', out);
    }

    format_default(opt, defaults, suffix, sizeof suffix);
    fprintf(out, "%s\n", suffix);
}

/* Groups are collected by scanning rather than by assuming the table is sorted
 * by group, so a row added anywhere lands under the right heading. Hidden rows
 * do not count, so a group starts at its first visible option and a group with
 * none at all never prints a heading. */
static bool group_seen_before(const cli_spec *spec, size_t index)
{
    for (size_t i = 0; i < index; i++)
        if (!spec->options[i].hidden &&
            strcmp(spec->options[i].group, spec->options[index].group) == 0)
            return true;

    return false;
}

/* A placeholder, with an ellipsis where it stands for any number of arguments. */
static const char *positional_form(const cli_positional *pos, char *out, size_t len)
{
    snprintf(out, len, "%s%s", pos->metavar, pos->variadic ? "..." : "");

    return out;
}

static void print_usage_line(const cli_spec *spec, FILE *out)
{
    fprintf(out, "usage: %s [options]", spec->program);

    for (size_t i = 0; i < spec->n_options; i++) {
        if (!spec->options[i].required)
            continue;

        if (spec->options[i].key)
            fprintf(out, " -%c %s", spec->options[i].key, spec->options[i].metavar);
        else
            fprintf(out, " --%s %s", spec->options[i].name, spec->options[i].metavar);
    }

    for (size_t i = 0; i < spec->n_positionals; i++) {
        char form[METAVAR_MAX];

        fprintf(out, spec->positionals[i].required ? " %s" : " [%s]",
                positional_form(&spec->positionals[i], form, sizeof form));
    }

    fputc('\n', out);
}

static void print_positionals(const cli_spec *spec, FILE *out)
{
    if (spec->n_positionals == 0)
        return;

    fprintf(out, "\nArguments\n");

    for (size_t i = 0; i < spec->n_positionals; i++) {
        char form[METAVAR_MAX];

        fprintf(out, "  %-*s %s\n", HELP_COLUMN,
                positional_form(&spec->positionals[i], form, sizeof form),
                spec->positionals[i].help);
    }
}

void cli_usage(const cli_spec *spec, FILE *out)
{
    fprintf(out, "%s %s -- %s\n\n", spec->program, spec->version, spec->summary);
    print_usage_line(spec, out);
    print_positionals(spec, out);

    for (size_t i = 0; i < spec->n_options; i++) {
        if (spec->options[i].hidden || group_seen_before(spec, i))
            continue;

        fprintf(out, "\n%s\n", spec->options[i].group);

        for (size_t j = i; j < spec->n_options; j++)
            if (!spec->options[j].hidden &&
                strcmp(spec->options[j].group, spec->options[i].group) == 0)
                print_option(out, &spec->options[j], spec->defaults);
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

static const char *type_name(cli_type type)
{
    switch (type) {
        case OPT_FLAG:   return "flag";
        case OPT_STRING: return "string";
        case OPT_SIZE:   return "size";
        case OPT_INT:    return "int";
        case OPT_DOUBLE: return "double";
        case OPT_ENUM:   return "enum";
    }

    return "unknown";
}

static void print_json_choices(FILE *out, const cli_option *opt)
{
    if (!opt->choices) {
        fputs("null", out);
        return;
    }

    fputc('[', out);
    for (const cli_choice *choice = opt->choices; choice->name; choice++) {
        if (choice != opt->choices)
            fputs(", ", out);
        print_json_string(out, choice->name);
    }
    fputc(']', out);
}

static void print_json_default(FILE *out, const cli_option *opt, const void *defaults)
{
    const char *field = (const char *)defaults + opt->offset;

    if (!stores_a_value(opt)) {
        fputs("null", out);
        return;
    }

    switch (opt->type) {
        case OPT_FLAG:   fputs(*(const bool *)field ? "true" : "false", out); break;
        case OPT_STRING: print_json_string(out, *(const char *const *)field); break;
        case OPT_SIZE:   fprintf(out, "%zu", *(const size_t *)field);         break;
        case OPT_INT:    fprintf(out, "%d", *(const int *)field);             break;
        case OPT_DOUBLE: fprintf(out, "%g", *(const double *)field);          break;
        case OPT_ENUM:
            for (const cli_choice *choice = opt->choices; choice->name; choice++)
                if (choice->value == *(const int *)field)
                    print_json_string(out, choice->name);
            break;
    }
}

static void print_json_bounds(FILE *out, const cli_option *opt)
{
    if (opt->type != OPT_SIZE && opt->type != OPT_INT && opt->type != OPT_DOUBLE) {
        fputs("      \"minimum\": null, \"maximum\": null\n", out);
        return;
    }

    fprintf(out, "      \"minimum\": %ld, \"maximum\": ", opt->minimum);

    if (opt->maximum == CLI_UNBOUNDED)
        fputs("null\n", out);
    else
        fprintf(out, "%ld\n", opt->maximum);
}

static void print_json_option(FILE *out, const cli_option *opt, const void *defaults,
                              bool last)
{
    char key[2] = { opt->key, '\0' };

    fputs("    {\n", out);
    fputs("      \"name\": ", out);           print_json_string(out, opt->name);
    fputs(",\n      \"short\": ", out);       print_json_string(out, opt->key ? key : NULL);
    fputs(",\n      \"group\": ", out);       print_json_string(out, opt->group);
    fprintf(out, ",\n      \"type\": \"%s\"", type_name(opt->type));
    fputs(",\n      \"metavar\": ", out);     print_json_string(out, opt->metavar);
    fputs(",\n      \"help\": ", out);        print_json_string(out, opt->help);
    fputs(",\n      \"detail\": ", out);      print_json_string(out, opt->detail);
    fprintf(out, ",\n      \"required\": %s", opt->required ? "true" : "false");
    fprintf(out, ",\n      \"hidden\": %s", opt->hidden ? "true" : "false");
    fputs(",\n      \"unset_label\": ", out); print_json_string(out, opt->unset_label);
    fputs(",\n      \"choices\": ", out);     print_json_choices(out, opt);
    fputs(",\n      \"default\": ", out);     print_json_default(out, opt, defaults);
    fputs(",\n", out);                        print_json_bounds(out, opt);
    fprintf(out, "    }%s\n", last ? "" : ",");
}

static void print_json_positional(FILE *out, const cli_positional *pos, bool last)
{
    fputs("    {\n", out);
    fputs("      \"name\": ", out);       print_json_string(out, pos->name);
    fputs(",\n      \"metavar\": ", out); print_json_string(out, pos->metavar);
    fputs(",\n      \"help\": ", out);    print_json_string(out, pos->help);
    fputs(",\n      \"detail\": ", out);  print_json_string(out, pos->detail);
    fprintf(out, ",\n      \"required\": %s", pos->required ? "true" : "false");
    fprintf(out, ",\n      \"variadic\": %s\n", pos->variadic ? "true" : "false");
    fprintf(out, "    }%s\n", last ? "" : ",");
}

void cli_dump_options(const cli_spec *spec, FILE *out)
{
    fputs("{\n", out);
    fputs("  \"program\": ", out); print_json_string(out, spec->program);
    fputs(",\n  \"version\": ", out); print_json_string(out, spec->version);
    fputs(",\n  \"summary\": ", out); print_json_string(out, spec->summary);
    fputs(",\n  \"options\": [\n", out);

    for (size_t i = 0; i < spec->n_options; i++)
        print_json_option(out, &spec->options[i], spec->defaults,
                          i + 1 == spec->n_options);

    fputs("  ],\n  \"positionals\": [\n", out);
    for (size_t i = 0; i < spec->n_positionals; i++)
        print_json_positional(out, &spec->positionals[i], i + 1 == spec->n_positionals);

    fputs("  ]\n}\n", out);
}

/* ------------------------------------------------------------------------ */
/* Parsing                                                                   */
/* ------------------------------------------------------------------------ */

static cli_status check_required_options(const cli_spec *spec, const bool *seen)
{
    for (size_t i = 0; i < spec->n_options; i++) {
        if (!spec->options[i].required || seen[i])
            continue;

        fprintf(stderr, "%s: missing required option --%s (%s)\n",
                spec->program, spec->options[i].name, spec->options[i].help);
        return CLI_ERROR;
    }

    return CLI_OK;
}

/* A variadic positional takes any number of arguments, but demands one where it
 * is required. */
static int fewest_arguments(const cli_spec *spec)
{
    int fewest = 0;

    for (size_t i = 0; i < spec->n_positionals; i++)
        if (!spec->positionals[i].variadic || spec->positionals[i].required)
            fewest++;

    return fewest;
}

static bool takes_any_number(const cli_spec *spec)
{
    return spec->n_positionals > 0 && spec->positionals[spec->n_positionals - 1].variadic;
}

static void report_positionals(const cli_spec *spec, int given)
{
    fprintf(stderr, "%s: expected ", spec->program);

    for (size_t i = 0; i < spec->n_positionals; i++) {
        char form[METAVAR_MAX];

        fprintf(stderr, "%s%s", i ? " " : "",
                positional_form(&spec->positionals[i], form, sizeof form));
    }

    fprintf(stderr, ", got %d argument%s\n", given, given == 1 ? "" : "s");
}

static void store_positional(const cli_positional *pos, void *args, char **argv,
                             int at, int argc)
{
    if (!pos->variadic) {
        *(const char **)((char *)args + pos->offset) = argv[at];
        return;
    }

    *(const char *const **)((char *)args + pos->offset)  = (const char *const *)(argv + at);
    *(size_t *)((char *)args + pos->count_offset)        = (size_t)(argc - at);
}

static cli_status take_positionals(const cli_spec *spec, int argc, char **argv,
                                   void *args)
{
    int given = argc - optind;
    int at    = optind;

    if (given < fewest_arguments(spec) ||
        (given > (int)spec->n_positionals && !takes_any_number(spec))) {
        report_positionals(spec, given);
        return CLI_ERROR;
    }

    for (size_t i = 0; i < spec->n_positionals; i++) {
        store_positional(&spec->positionals[i], args, argv, at, argc);
        at = spec->positionals[i].variadic ? argc : at + 1;
    }

    return CLI_OK;
}

/* Requests that the tables answer in full, before any argument is required to
 * be present. */
static bool answer(const cli_spec *spec, cli_action action)
{
    switch (action) {
        case CLI_SHOW_HELP:    cli_usage(spec, stdout);        return true;
        case CLI_SHOW_VERSION: printf("%s %s\n", spec->program, spec->version);
                               return true;
        case CLI_DUMP_OPTIONS: cli_dump_options(spec, stdout); return true;
        case CLI_STORE:        return false;
    }

    return false;
}

cli_status cli_parse(const cli_spec *spec, int argc, char **argv, void *args)
{
    size_t         shortopts_size = 2 * spec->n_options + 3;
    struct option *longopts  = calloc(spec->n_options + 1, sizeof *longopts);
    char          *shortopts = calloc(shortopts_size, sizeof *shortopts);
    bool          *seen      = calloc(spec->n_options, sizeof *seen);
    cli_action     requested = CLI_STORE;
    cli_status     status    = CLI_ERROR;
    int            found;

    if (!longopts || !shortopts || !seen) {
        fprintf(stderr, "%s: out of memory\n", spec->program);
        goto done;
    }

    memcpy(args, spec->defaults, spec->args_size);
    build_short_options(spec, shortopts, shortopts_size);
    build_long_options(spec, longopts);

    opterr = 0;
    optind = 1;

    /* getopt keeps its state in globals, which is what makes it unsafe to call
     * from more than one thread. Arguments are parsed before any thread exists,
     * and optind is reset above so that a second parse starts afresh.
     * NOLINTNEXTLINE(concurrency-mt-unsafe) */
    while ((found = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1) {
        const cli_option *opt = found < OPTION_ID_BASE
                              ? option_by_key(spec, (char)found)
                              : option_by_id(spec, found);

        if (found == MISSING_VALUE || !opt) {
            report_rejected_option(spec, argc, argv, found == MISSING_VALUE);
            goto done;
        }

        seen[opt - spec->options] = true;

        if (!stores_a_value(opt)) {
            requested = opt->action;
            continue;
        }

        if (assign(opt, args, optarg, spec->program) < 0)
            goto done;
    }

    if (answer(spec, requested)) {
        status = CLI_DONE;
        goto done;
    }

    if (check_required_options(spec, seen) != CLI_OK)
        goto done;

    status = take_positionals(spec, argc, argv, args);

done:
    free(seen);
    free(shortopts);
    free(longopts);
    return status;
}
