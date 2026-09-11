/* options.c -- cmuts csv's command line, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "options.h"

#include "manifest.h"

#include <stddef.h>
#include <stdio.h>

static void dump_layout(FILE *out)
{
    fmt_dump_reads(out, "cmuts csv", CMUTS_CSV_READS, CMUTS_CSV_N_READS);
}

static const cli_option OPTIONS[] = {
    {
        .group       = "Input",
        .name        = "fasta",
        .key         = 'f',
        .type        = OPT_STRING,
        .offset      = offsetof(csv_args, csv.fasta_path),
        .metavar     = "FASTA",
        .help        = "name the references from this file, in the order of the rows",
        .unset_label = "the row number",
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
        .help   = "describe the input format as JSON and exit",
        .hidden = true,
        .action = CLI_PRINT,
        .print  = dump_layout,
    },
};

static const cli_positional POSITIONALS[] = {
    {
        .name     = "input",
        .metavar  = "HDF5",
        .help     = "the cmuts output to convert",
        .offset   = offsetof(csv_args, csv.input_path),
        .required = true,
    },
};

csv_args csv_defaults(void)
{
    return (csv_args){ 0 };
}

cli_spec csv_spec(const csv_args *defaults)
{
    return (cli_spec){
        .program       = "cmuts csv",
        .version       = CMUTS_VERSION,
        .summary       = "write a cmuts output as comma separated values.",
        .options       = OPTIONS,
        .n_options     = sizeof OPTIONS / sizeof *OPTIONS,
        .positionals   = POSITIONALS,
        .n_positionals = sizeof POSITIONALS / sizeof *POSITIONALS,
        .defaults      = defaults,
        .args_size     = sizeof *defaults,
    };
}
