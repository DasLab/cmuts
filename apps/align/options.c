/* options.c -- cmuts align's command line, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "options.h"

#include <stddef.h>

#include "preset.h"

/* The default number of threads for merging, alignment and sorting. */
#define THREADS_DEFAULT 1

/* The memory one sorting thread may hold, in mebibytes. samtools sort writes a temporary
 * file every time a thread fills, so this is the default samtools itself uses. */
#define MEMORY_DEFAULT 768

static const cli_option OPTIONS[] = {
    {
        .group    = "Input and output",
        .name     = "fasta",
        .key      = 'f',
        .type     = OPT_STRING,
        .offset   = offsetof(align_args, mapping.fasta_path),
        .metavar  = "FASTA",
        .help     = "reference sequences",
        .required = true,
    },
    {
        .group    = "Input and output",
        .name     = "output",
        .key      = 'o',
        .type     = OPT_STRING,
        .offset   = offsetof(align_args, mapping.output_path),
        .metavar  = "BAM",
        .help     = "write sorted alignments to this file",
        .required = true,
    },
    {
        .group  = "Input and output",
        .name   = "overwrite",
        .type   = OPT_FLAG,
        .offset = offsetof(align_args, mapping.overwrite),
        .help   = "replace the output file if it already exists",
    },

    {
        .group    = "Alignment",
        .name     = "preset",
        .key      = 'x',
        .type     = OPT_ENUM,
        .offset   = offsetof(align_args, preset),
        .metavar  = "PRESET",
        .help     = "minimap2 preset for the sequencing technology",
        .label    = "Alignment Preset",
        .required = true,
        .choices  = PRESET_CHOICES,
    },

    {
        .group   = "Performance",
        .name    = "threads",
        .key     = 't',
        .type    = OPT_INT,
        .offset  = offsetof(align_args, mapping.threads),
        .metavar = "N",
        .help    = "threads for merging, alignment and sorting",
        .minimum = 1,
        .maximum = CLI_UNBOUNDED,
    },

    {
        .group   = "Performance",
        .name    = "memory",
        .key     = 'm',
        .type    = OPT_SIZE,
        .offset  = offsetof(align_args, mapping.memory),
        .metavar = "MiB",
        .help    = "memory each sorting thread holds before it writes to disk",
        .minimum = 1,
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
        .name     = "reads",
        .metavar  = "READS",
        .help     = "reads to align, as FASTQ or unaligned BAM",
        .offset   = offsetof(align_args, mapping.reads_paths[0]),
        .required = true,
    },
    {
        .name    = "mate",
        .metavar = "MATE",
        .help    = "the second file of a pair, merged with the first",
        .offset  = offsetof(align_args, mapping.reads_paths[1]),
    },
};

align_args align_defaults(void)
{
    return (align_args){
        .mapping = {
            .threads = THREADS_DEFAULT,
            .memory  = MEMORY_DEFAULT,
        },
        .preset  = PRESET_UNSET,
    };
}

/* Copies the preset, which the parser writes as an int, into the module's configuration. */
void align_take_arguments(align_args *args)
{
    args->mapping.preset = (preset_id)args->preset;
}

cli_spec align_spec(const align_args *defaults)
{
    return (cli_spec){
        .program       = "cmuts align",
        .version       = CMUTS_VERSION,
        .summary       = "align reads to a reference and sort the resulting alignments.",
        .options       = OPTIONS,
        .n_options     = sizeof OPTIONS / sizeof *OPTIONS,
        .positionals   = POSITIONALS,
        .n_positionals = sizeof POSITIONALS / sizeof *POSITIONALS,
        .defaults      = defaults,
        .args_size     = sizeof *defaults,
    };
}
