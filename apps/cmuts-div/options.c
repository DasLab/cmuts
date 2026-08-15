/* options.c -- cmuts-div's command line, as one table.
 * Rows use designated initializers, so that a field added to cli_option takes its
 * default without being spelled out in every row.
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 * */

#include "options.h"

#include <stddef.h>

static const cli_option OPTIONS[] = {
    {
        .group    = "Input and output",
        .name     = "output",
        .key      = 'o',
        .type     = OPT_STRING,
        .offset   = offsetof(div_args, divide.output_path),
        .metavar  = "HDF5",
        .help     = "write results to this file",
        .required = true,
    },
    {
        .group  = "Input and output",
        .name   = "overwrite",
        .type   = OPT_FLAG,
        .offset = offsetof(div_args, divide.overwrite),
        .help   = "replace the output file if it already exists",
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
};

static const cli_positional POSITIONALS[] = {
    {
        .name     = "rates",
        .metavar  = "RATES",
        .help     = "the reactivities to normalize",
        .offset   = offsetof(div_args, divide.rates_path),
        .required = true,
    },
    {
        .name     = "control",
        .metavar  = "CONTROL",
        .help     = "the denatured control",
        .offset   = offsetof(div_args, divide.control_path),
        .required = true,
    },
};

div_args div_defaults(void)
{
    return (div_args){ 0 };
}

cli_spec div_spec(const div_args *defaults)
{
    return (cli_spec){
        .program       = "cmuts-div",
        .version       = CMUTS_VERSION,
        .summary       = "divide a cmuts output by a denatured control.",
        .options       = OPTIONS,
        .n_options     = sizeof OPTIONS / sizeof *OPTIONS,
        .positionals   = POSITIONALS,
        .n_positionals = sizeof POSITIONALS / sizeof *POSITIONALS,
        .defaults      = defaults,
        .args_size     = sizeof *defaults,
    };
}
