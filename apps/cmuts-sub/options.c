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
      .detail = "Written in the same layout as the two inputs, so that anything "
                "reading a cmuts output reads this one as well.",
      .required = true },
    { .group = "Input and output", .name = "overwrite", .type = OPT_FLAG,
      .offset = offsetof(sub_args, subtract.overwrite),
      .help = "replace the output file if it already exists",
      .detail = "Without this a run stops rather than destroying a result that "
                "may have cost a great deal more to produce than the one being "
                "started." },

    { .group = "Subtraction", .name = "denatured", .key = 'd', .type = OPT_STRING,
      .offset = offsetof(sub_args, subtract.denatured_path), .metavar = "HDF5",
      .help = "normalize against a denatured control",
      .detail = "A sample denatured before it was treated, in which the reagent "
                "reaches every position. Dividing by its rates corrects each "
                "position for how readily the reagent modifies it at all. The "
                "result is then a ratio of two rates rather than a rate, and a "
                "position the control measured at zero is reported as unmeasured "
                "rather than as infinitely reactive.",
      .unset_label = "none" },

    { .group = "Subtraction", .name = "clip", .type = OPT_FLAG,
      .offset = offsetof(sub_args, subtract.clip),
      .help = "raise a negative reactivity to zero",
      .detail = "A rate below its background comes out negative, which no "
                "structure produces. Clipping reports those positions as "
                "unmodified instead, and discards how far below the background "
                "the treated sample fell. A position neither run measured stays "
                "unmeasured either way." },

    { .group = "Information", .name = "help", .key = 'h', .type = OPT_FLAG,
      .help = "show this help and exit", .action = CLI_SHOW_HELP },
    { .group = "Information", .name = "version", .key = 'V', .type = OPT_FLAG,
      .help = "show the version and exit", .action = CLI_SHOW_VERSION },
    { .group = "Information", .name = "dump-options", .type = OPT_FLAG,
      .help = "describe every argument as JSON and exit",
      .detail = "Intended for generating documentation and shell completions from "
                "the binary rather than from a separate description of it.",
      .hidden = true, .action = CLI_DUMP_OPTIONS },
};

#define SHARED_DETAIL \
    "Both inputs must have been counted against the same references, in the " \
    "same order. An output records no reference names -- the FASTA holds " \
    "those, and a row is identified by its position -- so what is checked is " \
    "that the two agree in the number of rows and in their width."

static const cli_positional POSITIONALS[] = {
    { .name = "treated", .metavar = "TREATED", .help = "the modified sample",
      .detail = "The run whose reactivities the background is taken off. "
                SHARED_DETAIL,
      .offset = offsetof(sub_args, subtract.treated_path), .required = true },
    { .name = "untreated", .metavar = "UNTREATED", .help = "the background",
      .detail = "The run standing for what is measured in the absence of the "
                "modifying reagent. " SHARED_DETAIL,
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
