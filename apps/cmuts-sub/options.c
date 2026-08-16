/* options.c -- cmuts-sub's command line, as one table.
 * Rows use designated initializers, so that a field added to cli_option takes its
 * default without being spelled out in every row.
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 * */

#include "options.h"

#include <stddef.h>

/* The datasets a run leaves behind, which are those every output holds: what it was
 * given, less the background it was given. */
static const out_written WRITTEN[] = {
    { .id = OUT_COVERAGE },
    { .id = OUT_REACTIVITY, .detail = "The mutation rate of the treated sample less that of the untreated one, so what remains is the signal the treatment added." },
    { .id = OUT_ERROR },
    { .id = OUT_LENGTHS },
    { .id = OUT_READS },
    { .id = OUT_REJECTED },
    { .id = OUT_UNMAPPED },
};

const out_manifest CMUTS_SUB_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };

static void dump_layout(FILE *out)
{
    out_dump_layout(out, "cmuts-sub", &CMUTS_SUB_WRITES);
}

static const cli_option OPTIONS[] = {
    {
        .group    = "Input and output",
        .name     = "output",
        .key      = 'o',
        .type     = OPT_STRING,
        .offset   = offsetof(sub_args, subtract.output_path),
        .metavar  = "HDF5",
        .help     = "write results to this file",
        .required = true,
    },
    {
        .group  = "Input and output",
        .name   = "overwrite",
        .type   = OPT_FLAG,
        .offset = offsetof(sub_args, subtract.overwrite),
        .help   = "replace the output file if it already exists",
    },

    {
        .group  = "Subtraction",
        .name   = "clip",
        .type   = OPT_FLAG,
        .offset = offsetof(sub_args, subtract.clip),
        .help   = "raise a negative reactivity to zero",
    },

    {
        .group  = "Information",
        .name   = "help",
        .key    = 'h',
        .type   = OPT_FLAG,
        .help   = "show this help and exit",
        .action = CLI_SHOW_HELP,
    },
    {
        .group  = "Information",
        .name   = "version",
        .key    = 'V',
        .type   = OPT_FLAG,
        .help   = "show the version and exit",
        .action = CLI_SHOW_VERSION,
    },
    {
        .group  = "Information",
        .name   = "dump-options",
        .type   = OPT_FLAG,
        .help   = "describe every argument as JSON and exit",
        .hidden = true,
        .action = CLI_DUMP_OPTIONS,
    },
    {
        .group  = "Information",
        .name   = "dump-layout",
        .type   = OPT_FLAG,
        .help   = "describe the output format as JSON and exit",
        .hidden = true,
        .action = CLI_PRINT,
        .print  = dump_layout,
    },
};

static const cli_positional POSITIONALS[] = {
    {
        .name     = "treated",
        .metavar  = "TREATED",
        .help     = "the modified sample",
        .offset   = offsetof(sub_args, subtract.treated_path),
        .required = true,
    },
    {
        .name     = "untreated",
        .metavar  = "UNTREATED",
        .help     = "the background",
        .offset   = offsetof(sub_args, subtract.untreated_path),
        .required = true,
    },
};

sub_args sub_defaults(void)
{
    return (sub_args){ 0 };
}

cli_spec sub_spec(const sub_args *defaults)
{
    return (cli_spec){
        .program       = "cmuts-sub",
        .version       = CMUTS_VERSION,
        .summary       = "subtract an untreated background from a cmuts output.",
        .options       = OPTIONS,
        .n_options     = sizeof OPTIONS / sizeof *OPTIONS,
        .positionals   = POSITIONALS,
        .n_positionals = sizeof POSITIONALS / sizeof *POSITIONALS,
        .defaults      = defaults,
        .args_size     = sizeof *defaults,
    };
}
