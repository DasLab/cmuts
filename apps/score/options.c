/* options.c -- cmuts score's command line, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "options.h"

#include <stddef.h>

static const cli_option OPTIONS[] = {
    {
        .group    = "Input",
        .name     = "fasta",
        .key      = 'f',
        .type     = OPT_STRING,
        .offset   = offsetof(score_args, score.fasta_path),
        .metavar  = "FASTA",
        .help     = "the references the input was counted against, in the order of its rows",
        .required = true,
    },
    {
        .group   = "Input",
        .name     = "structures",
        .key      = 's',
        .type     = OPT_STRING,
        .offset   = offsetof(score_args, score.structures_path),
        .metavar  = "FILE",
        .help     = "dot bracket records, matched to the references by name",
        .required = true,
    },

    {
        .group       = "Scoring",
        .name        = "bases",
        .key         = 'b',
        .type        = OPT_STRING,
        .offset      = offsetof(score_args, score.bases),
        .metavar     = "BASES",
        .help        = "score only the bases the reagent modifies",
        .unset_label = "every base",
    },
    {
        .group   = "Scoring",
        .name    = "min-coverage",
        .type    = OPT_DOUBLE,
        .offset  = offsetof(score_args, score.min_coverage),
        .metavar = "D",
        .help    = "reads a position needs before it is scored",
        .minimum = 0,
        .maximum = CLI_UNBOUNDED,
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
        .name     = "input",
        .metavar  = "HDF5",
        .help     = "the output being scored",
        .offset   = offsetof(score_args, score.input_path),
        .required = true,
    },
};

score_args score_defaults(void)
{
    return (score_args){ 0 };
}

cli_spec score_spec(const score_args *defaults)
{
    return (cli_spec){
        .program       = "cmuts score",
        .version       = CMUTS_VERSION,
        .summary       = "measure a cmuts output against a known structure.",
        .options       = OPTIONS,
        .n_options     = sizeof OPTIONS / sizeof *OPTIONS,
        .positionals   = POSITIONALS,
        .n_positionals = sizeof POSITIONALS / sizeof *POSITIONALS,
        .defaults      = defaults,
        .args_size     = sizeof *defaults,
    };
}
