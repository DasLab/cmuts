/* options.c -- cmuts-sub's command line, as one table.
 *
 * Rows use designated initializers, so that a field added to cli_option takes its default
 * rather than having to be spelled out in every row.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "options.h"

#include <stddef.h>

static const cli_option OPTIONS[] = {
    { .group = "Input and output", .name = "output", .key = 'o', .type = OPT_STRING,
      .offset = offsetof(sub_args, subtract.output_path), .metavar = "HDF5",
      .help = "write results to this file",
      .required = true },
    { .group = "Input and output", .name = "overwrite", .type = OPT_FLAG,
      .offset = offsetof(sub_args, subtract.overwrite),
      .help = "replace the output file if it already exists" },

    { .group = "Subtraction", .name = "denatured", .key = 'd', .type = OPT_STRING,
      .offset = offsetof(sub_args, subtract.denatured_path), .metavar = "HDF5",
      .help = "normalize against a denatured control",
      .unset_label = "none" },

    { .group = "Subtraction", .name = "clip", .type = OPT_FLAG,
      .offset = offsetof(sub_args, subtract.clip),
      .help = "raise a negative reactivity to zero" },

    { .group = "Information", .name = "help", .key = 'h', .type = OPT_FLAG,
      .help = "show this help and exit", .action = CLI_SHOW_HELP },
    { .group = "Information", .name = "version", .key = 'V', .type = OPT_FLAG,
      .help = "show the version and exit", .action = CLI_SHOW_VERSION },
    { .group = "Information", .name = "dump-options", .type = OPT_FLAG,
      .help = "describe every argument as JSON and exit",
      .hidden = true, .action = CLI_DUMP_OPTIONS },
};

static const cli_positional POSITIONALS[] = {
    { .name = "treated", .metavar = "TREATED", .help = "the modified sample",
      .offset = offsetof(sub_args, subtract.treated_path), .required = true },
    { .name = "untreated", .metavar = "UNTREATED", .help = "the background",
      .offset = offsetof(sub_args, subtract.untreated_path), .required = true },
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
