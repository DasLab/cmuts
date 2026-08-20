/* options.c -- cmuts sub's command line, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "options.h"

#include <stddef.h>

/* The datasets a run leaves behind, which are those every output holds: what it was
 * given, less the background it was given. */
static const char SUMMED[] = "Summed over the inputs.";

static const out_written WRITTEN[] = {
    { .id = OUT_REACTIVITY,
      .note = "Treated less untreated.",
      .origin = OUT_REQUIRED },
    { .id = OUT_ERROR,
      .note = "Propagated from the inputs." },
    { .id = OUT_COVERAGE,
      .note = SUMMED,
      .origin = OUT_REQUIRED },
    { .id = OUT_SEQUENCE },
    { .id = OUT_LENGTHS,   .note = SUMMED },
    { .id = OUT_READS,     .note = SUMMED },
    { .id = OUT_REJECTED,  .note = SUMMED },
    { .id = OUT_UNMAPPED,  .note = SUMMED },
};


const out_manifest CMUTS_SUB_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };

static void dump_layout(FILE *out)
{
    out_dump_layout(out, "cmuts sub", &CMUTS_SUB_WRITES);
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
        .program       = "cmuts sub",
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
