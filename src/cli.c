/* cli.c -- command line parsing
 *
 * Every argument a program accepts is described once, in its cli_spec. The short
 * option string, the getopt_long table, assignment, bounds checking, the usage
 * line, the help text and the JSON description all derive from it, so adding an
 * argument means adding a row.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "cli.h"

#include <errno.h>
#include <float.h>
#include <getopt.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Long options without a short form need an identifier getopt_long can return;
 * anything above the character range will do. */
#define OPTION_ID_BASE 256

/* What getopt returns for an option given without the value it needs. */
#define MISSING_VALUE ':'

#define REJECTED_MAX     64  /* an option as it was written on the command line */
#define METAVAR_MAX      32  /* a placeholder, with any ellipsis */
#define INVOCATION_MAX   64  /* both forms of an option, and its placeholder */
#define CEILING_MAX      32  /* an option's largest value, written out */
#define SET_LIST_MAX     64  /* the choices a set holds, comma separated */

/* The note giving an option's default. Its longest form wraps a full set list,
 * so it is sized from that. */
#define DEFAULT_NOTE_MAX (SET_LIST_MAX + sizeof(" (default )"))

/* A leading colon, a letter and possibly a colon per option, and the terminator. */
#define SHORTOPTS_SIZE(options) (2 + (2 * (options)))

/* The column an option's help begins in, shared with the arguments list so the
 * two cannot step apart. A floor and not a fixed width, since a metavar may be
 * wider than it leaves room for. */
#define HELP_COLUMN 28

/* Reports an option getopt would not accept, quoting it as the caller wrote it.
 * getopt records a short option in optopt but leaves a long one as the identifier
 * it was assigned, so a long one is read back from the word getopt stopped on. */
static void report_rejected_option(const cli_spec *spec, int argc, char **argv,
                                   bool needs_value)
{
    char rejected[REJECTED_MAX];

    if (optopt > 0 && optopt < OPTION_ID_BASE) {
        snprintf(rejected, sizeof rejected, "-%c", optopt);
    } else if (optind >= 1 && optind <= argc) {
        snprintf(rejected, sizeof rejected, "%s", argv[optind - 1]);
    } else {
        snprintf(rejected, sizeof rejected, "an option");
    }

    if (needs_value) {
        fprintf(stderr, "%s: %s needs a value\n", spec->program, rejected);
    } else {
        fprintf(stderr, "%s: unrecognized option %s\n", spec->program, rejected);
    }
}

static bool takes_argument(const cli_option *opt)
{
    return opt->type != OPT_FLAG;
}

static bool stores_a_value(const cli_option *opt)
{
    return opt->action == CLI_STORE;
}

/* ------------------------------------------------------------------------ */
/* Table lookups                                                             */
/* ------------------------------------------------------------------------ */

static const cli_option *option_by_key(const cli_spec *spec, char key)
{
    for (size_t i = 0; i < spec->n_options; i++) {
        if (spec->options[i].key == key) {
            return &spec->options[i];
        }
    }

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

/* Builds the short option string. The leading colon has getopt report a missing
 * value separately from an unknown option, so they can get different messages. The guard
 * leaves room for an option and the terminator, which SHORTOPTS_SIZE provides. */
static void build_short_options(const cli_spec *spec, char *out, size_t size)
{
    size_t n = 0;

    out[n++] = ':';

    for (size_t i = 0; i < spec->n_options && n + 2 < size; i++) {
        if (!spec->options[i].key) {
            continue;
        }

        out[n++] = spec->options[i].key;
        if (takes_argument(&spec->options[i])) {
            out[n++] = ':';
        }
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
    for (const cli_choice *choice = opt->choices; choice->name; choice++) {
        fprintf(out, "%s%s", choice == opt->choices ? "" : "|", choice->name);
    }
}

/* Takes the choice named by the first len characters, so that a word within a longer
 * string is matched without a copy of it. */
static int parse_choice_n(const cli_option *opt, const char *text, size_t len,
                          const char *program, int *out)
{
    for (const cli_choice *choice = opt->choices; choice->name; choice++) {
        if (strncmp(choice->name, text, len) == 0 && choice->name[len] == '\0') {
            *out = choice->value;
            return 0;
        }
    }

    fprintf(stderr, "%s: --%s: \"%.*s\" is not one of ", program, opt->name, (int)len,
            text);
    print_choices(stderr, opt);
    fputc('\n', stderr);

    return -1;
}

static int parse_choice(const cli_option *opt, const char *text, const char *program,
                        int *out)
{
    return parse_choice_n(opt, text, strlen(text), program, out);
}

/* Parses every choice a comma-separated list names, OR'd together.
 *
 * The empty subset is a choice worth zero. Asking for it alongside another is a
 * contradiction and not a preference between the two, so it is refused. */
static int parse_set(const cli_option *opt, const char *text, const char *program,
                     int *out)
{
    int  chosen = 0;
    bool empty  = false;

    for (const char *token = text; token; ) {
        const char *comma = strchr(token, ',');
        size_t      len   = comma ? (size_t)(comma - token) : strlen(token);
        int         choice;

        if (parse_choice_n(opt, token, len, program, &choice) < 0) {
            return -1;
        }

        empty  |= choice == 0;
        chosen |= choice;
        token   = comma ? comma + 1 : NULL;
    }

    if (empty && chosen != 0) {
        fprintf(stderr, "%s: --%s: \"%s\" combines the empty choice with another\n",
                program, opt->name, text);
        return -1;
    }

    *out = chosen;
    return 0;
}

/* Returns the largest whole number an option will take. One with no ceiling of its own
 * is bounded by what its destination can hold, so no value is accepted that the field
 * would silently truncate. */
static long integer_ceiling(const cli_option *opt)
{
    if (opt->maximum != CLI_UNBOUNDED) {
        return opt->maximum;
    }

    return opt->type == OPT_INT ? INT32_MAX : LONG_MAX;
}

static double real_ceiling(const cli_option *opt)
{
    return opt->maximum != CLI_UNBOUNDED ? (double)opt->maximum : DBL_MAX;
}

/* Returns the largest value an option will take, written out. */
static const char *ceiling_text(const cli_option *opt, char *out, size_t size)
{
    if (opt->type == OPT_DOUBLE) {
        snprintf(out, size, "%g", real_ceiling(opt));
    } else {
        snprintf(out, size, "%ld", integer_ceiling(opt));
    }

    return out;
}

/* Reports a value outside an option's bounds, echoing the caller's text rather
 * than the number parsed from it. A range is quoted only where the row declared
 * one; an option with no ceiling of its own is reported against whichever single
 * bound it missed. */
static void report_out_of_range(const cli_option *opt, const char *text,
                                const char *program, bool below)
{
    char limit[CEILING_MAX];

    if (opt->maximum != CLI_UNBOUNDED) {
        fprintf(stderr, "%s: --%s: %s is outside %ld..%ld\n",
                program, opt->name, text, opt->minimum, opt->maximum);
    } else if (below) {
        fprintf(stderr, "%s: --%s: %s is below %ld\n",
                program, opt->name, text, opt->minimum);
    } else {
        fprintf(stderr, "%s: --%s: %s is above %s\n",
                program, opt->name, text, ceiling_text(opt, limit, sizeof limit));
    }
}

static int parse_number(const cli_option *opt, const char *text, const char *program,
                        long *out)
{
    char *end = NULL;
    long  n;

    errno = 0;
    n     = strtol(text, &end, 10);

    if (*text == '\0' || *end != '\0') {
        fprintf(stderr, "%s: --%s: \"%s\" is not a number\n", program, opt->name, text);
        return -1;
    }

    /* strtol saturates at LONG_MIN or LONG_MAX on a value too wide for a long, which
     * an unbounded ceiling would otherwise accept. The saturated value carries the
     * direction it overflowed in. */
    if (errno == ERANGE || n < opt->minimum || n > integer_ceiling(opt)) {
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

/* Adds one value to a repeatable option's array, refusing the value that would overrun
 * it. */
static int append(const cli_option *opt, void *args, const char *value,
                  const char *program)
{
    const char **into  = (const char **)((char *)args + opt->offset);
    size_t      *count = (size_t *)((char *)args + opt->count_offset);

    if (*count == opt->capacity) {
        fprintf(stderr, "%s: --%s: at most %zu may be given\n", program, opt->name,
                opt->capacity);
        return -1;
    }

    into[(*count)++] = value;
    return 0;
}

static int assign(const cli_option *opt, void *args, const char *value,
                  const char *program)
{
    char  *field = (char *)args + opt->offset;
    long   number;
    double real;
    int    choice;

    if (opt->repeatable) {
        return append(opt, args, value, program);
    }

    switch (opt->type) {
        case OPT_DOUBLE:
            if (parse_double(opt, value, program, &real) < 0) {
                return -1;
            }
            *(double *)field = real;
            return 0;

        case OPT_ENUM:
            if (parse_choice(opt, value, program, &choice) < 0) {
                return -1;
            }
            *(int *)field = choice;
            return 0;

        case OPT_SET:
            if (parse_set(opt, value, program, &choice) < 0) {
                return -1;
            }
            *(int *)field = choice;
            return 0;

        case OPT_FLAG:
            *(bool *)field = true;
            return 0;

        case OPT_STRING:
            *(const char **)field = value;
            return 0;

        case OPT_SIZE:
            if (parse_number(opt, value, program, &number) < 0) {
                return -1;
            }
            *(size_t *)field = (size_t)number;
            return 0;

        case OPT_INT:
            if (parse_number(opt, value, program, &number) < 0) {
                return -1;
            }
            *(int *)field = (int)number;
            return 0;
    }

    return -1;
}

/* ------------------------------------------------------------------------ */
/* Help                                                                      */
/* ------------------------------------------------------------------------ */

/* Returns the names of every choice a set holds, comma separated, which is the spelling
 * the option itself takes. */
static void format_set(const cli_option *opt, int value, char *out, size_t size)
{
    size_t used = 0;

    out[0] = '\0';

    for (const cli_choice *choice = opt->choices; choice->name; choice++) {
        bool held = choice->value ? (value & choice->value) == choice->value
                                  : value == 0;
        int  n;

        if (!held) {
            continue;
        }

        n = snprintf(out + used, size - used, "%s%s", used ? "," : "", choice->name);

        if (n < 0 || (size_t)n >= size - used) {
            return;
        }

        used += (size_t)n;
    }
}

/* Renders the note giving an option's default, read from the spec's defaults so
 * that the help cannot advertise a value the program does not use. */
static void format_default(const cli_option *opt, const void *defaults,
                           char *out, size_t size)
{
    const char *field = (const char *)defaults + opt->offset;

    out[0] = '\0';

    if (!stores_a_value(opt)) {
        return;
    }

    /* An optional setting shows its label, not the sentinel value that means
     * "not applied". */
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
            if (*(const char *const *)field) {
                snprintf(out, size, " (default %s)", *(const char *const *)field);
            }
            break;
        case OPT_ENUM:
            for (const cli_choice *choice = opt->choices; choice->name; choice++) {
                if (choice->value == *(const int *)field) {
                    snprintf(out, size, " (default %s)", choice->name);
                }
            }
            break;
        case OPT_SET: {
            char names[SET_LIST_MAX];

            format_set(opt, *(const int *)field, names, sizeof names);
            snprintf(out, size, " (default %s)", names);
            break;
        }
        default:
            break;
    }
}

/* Returns how an option is written on its help line: the short form where it has one,
 * the long form always, and its placeholder. Both the width pass and the print read it
 * from here, so neither can measure one thing and write another. */
static const char *option_form(const cli_option *opt, char *out, size_t len)
{
    int n;

    if (opt->key) {
        n = snprintf(out, len, "-%c, --%s", opt->key, opt->name);
    } else {
        n = snprintf(out, len, "    --%s", opt->name);
    }

    if (opt->metavar && n > 0 && (size_t)n < len) {
        snprintf(out + n, len - (size_t)n, " %s", opt->metavar);
    }

    return out;
}

static void print_option(FILE *out, const cli_option *opt, const void *defaults,
                         int column)
{
    char invocation[INVOCATION_MAX];
    char suffix[DEFAULT_NOTE_MAX];

    option_form(opt, invocation, sizeof invocation);

    fprintf(out, "  %-*s %s", column, invocation, opt->help);

    if (opt->choices) {
        fputs(" (", out);
        print_choices(out, opt);
        fputc(')', out);
    }

    format_default(opt, defaults, suffix, sizeof suffix);
    fprintf(out, "%s\n", suffix);
}

/* Returns whether an earlier visible option shares this one's group. Found by
 * scanning the table, which is not assumed sorted by group, so a row added
 * anywhere lands under the right heading and a group of only hidden rows prints
 * none. */
static bool group_seen_before(const cli_spec *spec, size_t index)
{
    for (size_t i = 0; i < index; i++) {
        if (!spec->options[i].hidden &&
            strcmp(spec->options[i].group, spec->options[index].group) == 0) {
            return true;
        }
    }

    return false;
}

/* Returns a placeholder, with an ellipsis where it stands for any number of arguments. */
static const char *positional_form(const cli_positional *pos, char *out, size_t len)
{
    snprintf(out, len, "%s%s", pos->metavar, pos->variadic ? "..." : "");

    return out;
}

static void print_usage_line(const cli_spec *spec, FILE *out)
{
    fprintf(out, "usage: %s [options]", spec->program);

    for (size_t i = 0; i < spec->n_options; i++) {
        if (!spec->options[i].required) {
            continue;
        }

        if (spec->options[i].key) {
            fprintf(out, " -%c %s", spec->options[i].key, spec->options[i].metavar);
        } else {
            fprintf(out, " --%s %s", spec->options[i].name, spec->options[i].metavar);
        }
    }

    for (size_t i = 0; i < spec->n_positionals; i++) {
        char form[METAVAR_MAX];

        fprintf(out, spec->positionals[i].required ? " %s" : " [%s]",
                positional_form(&spec->positionals[i], form, sizeof form));
    }

    fputc('\n', out);
}

/* Returns the widest invocation the help will print, so that no invocation runs into its own
 * description. Hidden options are left out and never printed. */
static int help_column(const cli_spec *spec)
{
    size_t widest = HELP_COLUMN;

    for (size_t i = 0; i < spec->n_options; i++) {
        char   form[INVOCATION_MAX];
        size_t len;

        if (spec->options[i].hidden) {
            continue;
        }

        len = strlen(option_form(&spec->options[i], form, sizeof form));
        widest = len > widest ? len : widest;
    }

    for (size_t i = 0; i < spec->n_positionals; i++) {
        char   form[METAVAR_MAX];
        size_t len = strlen(positional_form(&spec->positionals[i],
                                            form, sizeof form));

        widest = len > widest ? len : widest;
    }

    return (int)widest;
}

static void print_positionals(const cli_spec *spec, FILE *out, int column)
{
    if (spec->n_positionals == 0) {
        return;
    }

    fprintf(out, "\nArguments\n");

    for (size_t i = 0; i < spec->n_positionals; i++) {
        char form[METAVAR_MAX];

        fprintf(out, "  %-*s %s\n", column,
                positional_form(&spec->positionals[i], form, sizeof form),
                spec->positionals[i].help);
    }
}

void cli_usage(const cli_spec *spec, FILE *out)
{
    int column = help_column(spec);

    fprintf(out, "%s %s -- %s\n\n", spec->program, spec->version, spec->summary);
    print_usage_line(spec, out);
    print_positionals(spec, out, column);

    for (size_t i = 0; i < spec->n_options; i++) {
        if (spec->options[i].hidden || group_seen_before(spec, i)) {
            continue;
        }

        fprintf(out, "\n%s\n", spec->options[i].group);

        for (size_t j = i; j < spec->n_options; j++) {
            if (!spec->options[j].hidden &&
                strcmp(spec->options[j].group, spec->options[i].group) == 0) {
                print_option(out, &spec->options[j], spec->defaults, column);
            }
        }
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
                if (*p < 0x20) {
                    fprintf(out, "\\u%04x", *p);
                } else {
                    fputc(*p, out);
                }
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
        case OPT_SET:    return "set";
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
        if (choice != opt->choices) {
            fputs(", ", out);
        }
        print_json_string(out, choice->name);
    }
    fputc(']', out);
}

/* JSON has no spelling for a value that is not finite, so such a default is described as
 * having no default. */
static void print_json_double(FILE *out, double value)
{
    if (isfinite(value)) {
        fprintf(out, "%g", value);
    } else {
        fputs("null", out);
    }
}

static void print_json_default(FILE *out, const cli_option *opt, const void *defaults)
{
    const char *field = (const char *)defaults + opt->offset;

    if (!stores_a_value(opt) || opt->repeatable) {
        fputs("null", out);
        return;
    }

    switch (opt->type) {
        case OPT_FLAG:   fputs(*(const bool *)field ? "true" : "false", out); break;
        case OPT_STRING: print_json_string(out, *(const char *const *)field); break;
        case OPT_SIZE:   fprintf(out, "%zu", *(const size_t *)field);         break;
        case OPT_INT:    fprintf(out, "%d", *(const int *)field);             break;
        case OPT_DOUBLE: print_json_double(out, *(const double *)field);      break;
        case OPT_ENUM:
            for (const cli_choice *choice = opt->choices; choice->name; choice++) {
                if (choice->value == *(const int *)field) {
                    print_json_string(out, choice->name);
                }
            }
            break;
        case OPT_SET: {
            char names[SET_LIST_MAX];

            format_set(opt, *(const int *)field, names, sizeof names);
            print_json_string(out, names);
            break;
        }
    }
}

static void print_json_bounds(FILE *out, const cli_option *opt)
{
    if (opt->type != OPT_SIZE && opt->type != OPT_INT && opt->type != OPT_DOUBLE) {
        fputs("      \"minimum\": null, \"maximum\": null\n", out);
        return;
    }

    fprintf(out, "      \"minimum\": %ld, \"maximum\": ", opt->minimum);

    if (opt->maximum == CLI_UNBOUNDED) {
        fputs("null\n", out);
    } else {
        fprintf(out, "%ld\n", opt->maximum);
    }
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
    fprintf(out, ",\n      \"required\": %s", opt->required ? "true" : "false");
    fprintf(out, ",\n      \"hidden\": %s", opt->hidden ? "true" : "false");
    fprintf(out, ",\n      \"repeatable\": %s", opt->repeatable ? "true" : "false");
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

    for (size_t i = 0; i < spec->n_options; i++) {
        print_json_option(out, &spec->options[i], spec->defaults,
                          i + 1 == spec->n_options);
    }

    fputs("  ],\n  \"positionals\": [\n", out);
    for (size_t i = 0; i < spec->n_positionals; i++) {
        print_json_positional(out, &spec->positionals[i], i + 1 == spec->n_positionals);
    }

    fputs("  ]\n}\n", out);
}

/* ------------------------------------------------------------------------ */
/* Parsing                                                                   */
/* ------------------------------------------------------------------------ */

static cli_status check_required_options(const cli_spec *spec, const bool *seen)
{
    for (size_t i = 0; i < spec->n_options; i++) {
        if (!spec->options[i].required || seen[i]) {
            continue;
        }

        fprintf(stderr, "%s: missing required option --%s (%s)\n",
                spec->program, spec->options[i].name, spec->options[i].help);
        return CLI_ERROR;
    }

    return CLI_OK;
}

/* Returns how many positional arguments the spec demands. A variadic one takes any
 * number but demands one where it is required. */
static int fewest_arguments(const cli_spec *spec)
{
    int fewest = 0;

    for (size_t i = 0; i < spec->n_positionals; i++) {
        if (!spec->positionals[i].variadic || spec->positionals[i].required) {
            fewest++;
        }
    }

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

/* Handles a request the tables can satisfy without arguments, returning whether one
 * was made. Called before any argument is required to be present. */
static bool answer(const cli_spec *spec, const cli_option *opt)
{
    if (!opt) {
        return false;
    }

    switch (opt->action) {
        case CLI_SHOW_HELP:    cli_usage(spec, stdout);        return true;
        case CLI_SHOW_VERSION: printf("%s %s\n", spec->program, spec->version);
                               return true;
        case CLI_DUMP_OPTIONS: cli_dump_options(spec, stdout); return true;
        case CLI_PRINT:        opt->print(stdout);             return true;
        case CLI_STORE:        return false;
    }

    return false;
}

cli_status cli_parse(const cli_spec *spec, int argc, char **argv, void *args)
{
    size_t         shortopts_size = SHORTOPTS_SIZE(spec->n_options);
    struct option *longopts  = calloc(spec->n_options + 1, sizeof *longopts);
    char          *shortopts = calloc(shortopts_size, sizeof *shortopts);
    bool          *seen      = calloc(spec->n_options, sizeof *seen);
    const cli_option *requested = NULL;
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

    /* getopt keeps its state in globals. Arguments are parsed before any thread
     * exists, and optind is reset above so that a second parse starts afresh.
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
            requested = opt;
            continue;
        }

        if (assign(opt, args, optarg, spec->program) < 0) {
            goto done;
        }
    }

    if (answer(spec, requested)) {
        status = CLI_DONE;
        goto done;
    }

    if (check_required_options(spec, seen) != CLI_OK) {
        goto done;
    }

    status = take_positionals(spec, argc, argv, args);

done:
    free(seen);
    free(shortopts);
    free(longopts);
    return status;
}
