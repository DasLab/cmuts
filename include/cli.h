/* cli.h -- command line parsing, driven by a table of options.
 *
 * A program describes its command line once, as a cli_spec, and everything else derives
 * from it: parsing, bounds checking, the usage line, grouped help, and a JSON description
 * for generating documentation and shell completions. Nothing here knows what any
 * particular program's options mean.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

/* The maximum of an option with no ceiling of its own, bounded only by what its
 * destination can hold. It is also the widest number a command line can be read into, so a
 * row naming it outright means the same thing. */
#define CLI_UNBOUNDED LONG_MAX

/* The type of a row's destination. It must match the C type of the field at that offset,
 * the value being written through a pointer of exactly this type. */
typedef enum {
    OPT_FLAG,    /* takes no argument; sets a bool */
    OPT_STRING,
    OPT_SIZE,
    OPT_INT,
    OPT_DOUBLE,  /* bounds are still written as whole numbers */
    OPT_ENUM,    /* one of a named set of values; stores an int */
} cli_type;

/* One accepted value of an OPT_ENUM option. A choice list ends with a NULL name. */
typedef struct {
    const char *name;
    int         value;
} cli_choice;

/* What an option does besides storing a value. The three that answer and exit are declared
 * rather than recognized by name, so they need no field of their own and any program may
 * have them. */
typedef enum {
    CLI_STORE,
    CLI_SHOW_HELP,
    CLI_SHOW_VERSION,
    CLI_DUMP_OPTIONS,
} cli_action;

typedef struct {
    const char       *group;    /* heading this option appears under */
    const char       *name;     /* long form */
    char              key;      /* short form, or 0 for none */
    cli_type          type;
    size_t            offset;   /* destination within the args struct */
    const char       *metavar;  /* argument placeholder; NULL when it takes none */
    const char       *help;     /* one line, for the help output */
    const char       *detail;   /* paragraph for manual pages; NULL to reuse help */
    bool              required;
    /* The word the help prints in place of a default, for an option that need not be
     * applied at all. NULL where every value is a real setting. */
    const char       *unset_label;
    long              minimum;  /* bounds for the numeric types */
    long              maximum;  /* CLI_UNBOUNDED where only the floor binds */
    bool              hidden;   /* kept out of the help, still described by JSON */
    const cli_choice *choices;  /* accepted values, for OPT_ENUM */
    cli_action        action;
} cli_option;

typedef struct {
    const char *name;     /* identifier for documentation and completions */
    const char *metavar;  /* placeholder in the usage line */
    const char *help;
    const char *detail;
    size_t      offset;
    bool        required;
    /* Takes every remaining argument rather than one, so it must come last. The array goes
     * to offset and its length to count_offset. */
    bool        variadic;
    size_t      count_offset;
} cli_positional;

/* Everything one program's command line consists of. */
typedef struct {
    const char           *program;
    const char           *version;
    const char           *summary;  /* one line, shown above the usage */
    const cli_option     *options;
    size_t                n_options;
    const cli_positional *positionals;
    size_t                n_positionals;
    const void           *defaults;  /* an args struct with defaults filled in */
    size_t                args_size;
} cli_spec;

typedef enum {
    CLI_OK,     /* arguments parsed; carry on */
    CLI_DONE,   /* the request was answered in full; exit successfully */
    CLI_ERROR,  /* usage error, already reported */
} cli_status;

/* Fills args, which must be spec->args_size bytes, starting from the spec's defaults. */
cli_status cli_parse(const cli_spec *spec, int argc, char **argv, void *args);

void cli_usage(const cli_spec *spec, FILE *out);
void cli_dump_options(const cli_spec *spec, FILE *out);
