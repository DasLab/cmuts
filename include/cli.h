/* cli.h -- command line parsing
 *
 * A program describes its command line once, as a cli_spec, and everything else derives
 * from it: parsing, bounds checking, the usage line, grouped help, and a JSON description
 * for generating documentation and shell completions.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

/* The maximum of an option with no ceiling of its own, bounded only by what its
 * destination can hold. */
#define CLI_UNBOUNDED LONG_MAX

/* An OPT_SET option that accepts the empty set takes this name for it. */
#define CLI_NONE "none"

/* The type of a row's destination. It must match the C type of the field at that offset,
 * the value being written through a pointer of exactly this type. */
typedef enum {
    OPT_FLAG,    /* takes no argument; sets a bool */
    OPT_STRING,
    OPT_SIZE,
    OPT_INT,
    OPT_DOUBLE,  /* bounds are still written as whole numbers */
    OPT_ENUM,    /* one of a named set of values; stores an int */
    /* An OPT_SET option takes one or more values from a named set, and stores them OR'd
     * together. Each choice must therefore have a nonzero value with a bit of its own.
     * An option that sets accepts_none also takes CLI_NONE, which stores zero. CLI_NONE
     * cannot be combined with any other choice. */
    OPT_SET,
} cli_type;

/* One accepted value of an OPT_ENUM or OPT_SET option. A choice list ends with a NULL
 * name. The label is a human-readable name for generated interfaces; NULL where the
 * choice name reads well enough on its own. */
typedef struct {
    const char *name;
    int         value;
    const char *label;
} cli_choice;

/* What an option does besides storing a value. Those that print and exit need no field of
 * their own. CLI_PRINT calls the row's print function. */
typedef enum {
    CLI_STORE,
    CLI_SHOW_HELP,
    CLI_SHOW_VERSION,
    CLI_DUMP_OPTIONS,
    CLI_PRINT,
} cli_action;

/* One option applies only when another option is set to one of these choices. An option
 * given outside them is refused, so a setting that does nothing is reported. The option
 * named must be an OPT_ENUM. */
typedef struct {
    const char *option;   /* long name of the option that governs this one */
    const char *choices;  /* the choices it must be set to, comma separated */
} cli_condition;

typedef struct {
    const char       *group;    /* heading this option appears under */
    const char       *name;     /* long form */
    char              key;      /* short form, or 0 for none */
    cli_type          type;
    size_t            offset;   /* destination within the args struct */
    const char       *metavar;  /* argument placeholder; NULL when it takes none */
    const char       *help;     /* one line, for the help output */
    /* Human-readable name for generated interfaces such as the web form. NULL where the
     * option name reads well enough on its own. */
    const char       *label;
    bool              required;
    /* The word the help prints in place of a default, for an option that need not be
     * applied at all. NULL where every value is a real setting. */
    const char       *unset_label;
    long              minimum;  /* bounds for the numeric types */
    long              maximum;  /* CLI_UNBOUNDED where only the floor binds */
    bool              hidden;   /* kept out of the help, still described by JSON */
    /* Takes a value every time it appears, rather than the last one winning. The values
     * go to offset as an array of const char *, their count to count_offset, and capacity
     * is how many that array holds. Only OPT_STRING may repeat. */
    bool              repeatable;
    size_t            count_offset;
    size_t            capacity;
    const cli_choice *choices;  /* accepted values, for OPT_ENUM and OPT_SET */
    bool              accepts_none;  /* An OPT_SET option takes CLI_NONE if this is true. */
    cli_condition     applies_when;  /* the condition this option applies under; empty
                                        where it applies always */
    cli_action        action;
    void            (*print)(FILE *out);  /* what CLI_PRINT calls */
} cli_option;

typedef struct {
    const char *name;     /* identifier for documentation and completions */
    const char *metavar;  /* placeholder in the usage line */
    const char *help;
    size_t      offset;
    bool        required;
    /* Takes every remaining argument, so it must come last. The array goes to offset and
     * its length to count_offset. */
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
