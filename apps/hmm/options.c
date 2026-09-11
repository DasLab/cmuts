/* options.c -- cmuts hmm's command line, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "options.h"

#include "manifest.h"

#include <stddef.h>

#include "filter.h"
#include "format.h"
#include "pairs.h"
#include "params.h"
#include "phred.h"

static const cli_choice VERIFY_CHOICES[] = {
    { "name",     REFSEQ_VERIFY_NAME     },
    { "checksum", REFSEQ_VERIFY_CHECKSUM },
    { "none",     0                      },
    { NULL,       0                      },
};

static const cli_choice PAIRWISE_CHOICES[] = {
    { "correlation", PAIRS_CORRELATION },
    { "conditional", PAIRS_CONDITIONAL },
    { "none",        0                 },
    { NULL,          0                 },
};

static const cli_choice STRAND_CHOICES[] = {
    { "forward", FILTER_STRAND_FORWARD },
    { "reverse", FILTER_STRAND_REVERSE },
    { NULL,      0                     },
};

static void dump_layout(FILE *out)
{
    fmt_dump_layout(out, "cmuts hmm", &CMUTS_HMM_WRITES);
}

static const cli_option OPTIONS[] = {
    {
        .group    = "Input and output",
        .name     = "fasta",
        .key      = 'f',
        .type     = OPT_STRING,
        .offset   = offsetof(cli_args, pipeline.fasta_path),
        .metavar  = "FASTA",
        .help     = "reference sequences",
        .required = true,
    },
    {
        .group    = "Input and output",
        .name     = "output",
        .key      = 'o',
        .type     = OPT_STRING,
        .offset   = offsetof(cli_args, pipeline.output_path),
        .metavar  = "HDF5",
        .help     = "write results to this file",
        .required = true,
    },
    {
        .group  = "Input and output",
        .name   = "overwrite",
        .type   = OPT_FLAG,
        .offset = offsetof(cli_args, pipeline.overwrite),
        .help   = "replace the output file if it already exists",
    },
    {
        .group   = "Input and output",
        .name    = "verify",
        .type    = OPT_SET,
        .offset  = offsetof(cli_args, pipeline.verify),
        .metavar = "CHECKS",
        .help    = "identity checks to make against the FASTA",
        .choices = VERIFY_CHOICES,
    },

    {
        .group   = "Filtering",
        .name    = "min-mapq",
        .key     = 'q',
        .type    = OPT_INT,
        .offset  = offsetof(cli_args, pipeline.filter_config.min_mapq),
        .metavar = "N",
        .help    = "reject alignments below this mapping quality",
        .label   = "Minimum Mapping Quality",
        .minimum = 0,
        .maximum = FILTER_MAPQ_UNAVAILABLE - 1,
    },
    {
        .group       = "Filtering",
        .name        = "min-length",
        .type        = OPT_INT,
        .offset      = offsetof(cli_args, pipeline.filter_config.min_length),
        .metavar     = "N",
        .help        = "reject reads shorter than this",
        .label       = "Minimum Read Length",
        .unset_label = "no limit",
        .minimum     = 0,
        .maximum     = CLI_UNBOUNDED,
    },
    {
        .group       = "Filtering",
        .name        = "max-length",
        .type        = OPT_INT,
        .offset      = offsetof(cli_args, pipeline.filter_config.max_length),
        .metavar     = "N",
        .help        = "reject reads longer than this",
        .label       = "Maximum Read Length",
        .unset_label = "no limit",
        .minimum     = 0,
        .maximum     = CLI_UNBOUNDED,
    },
    {
        .group  = "Filtering",
        .name   = "drop-supplementary",
        .type   = OPT_FLAG,
        .offset = offsetof(cli_args, pipeline.filter_config.drop_supplementary),
        .help   = "reject all but the primary piece of a split read",
        .label  = "Drop Supplementary Alignments",
    },
    {
        .group   = "Filtering",
        .name    = "strand",
        .key     = 's',
        .type    = OPT_SET,
        .offset  = offsetof(cli_args, pipeline.filter_config.strand),
        .metavar = "STRANDS",
        .help    = "reject alignments not on these strands",
        .label   = "Strands",
        .choices = STRAND_CHOICES,
    },

    {
        .group   = "Counting",
        .name    = "band",
        .type    = OPT_INT,
        .offset  = offsetof(cli_args, pipeline.tally_config.band),
        .metavar = "N",
        .help    = "reference positions the pair HMM may look either side of the CIGAR",
        .label   = "HMM Bandwidth",
        .minimum = 0,
        .maximum = CLI_UNBOUNDED,
    },
    {
        .group   = "Counting",
        .name    = "min-phred",
        .type    = OPT_INT,
        .offset  = offsetof(cli_args, pipeline.tally_config.min_phred),
        .metavar = "Q",
        .help    = "assign the maximum sequencing error to bases below this quality",
        .label   = "Minimum Base Quality",
        .minimum = 0,
        .maximum = PHRED_MAX,
    },
    {
        .group   = "Counting",
        .name    = "pairwise",
        .type    = OPT_SET,
        .offset  = offsetof(cli_args, pipeline.pairwise),
        .metavar = "STATS",
        .help    = "write these statistics of how often two positions are modified together",
        .label   = "Pairwise Statistics",
        .choices = PAIRWISE_CHOICES,
    },
    {
        .group   = "Counting",
        .name    = "min-depth",
        .type    = OPT_DOUBLE,
        .offset  = offsetof(cli_args, pipeline.rate_config.min_depth),
        .metavar = "D",
        .help    = "write NaN rates and errors for positions below this depth",
        .label   = "Minimum Depth",
        .minimum = 0,
        .maximum = CLI_UNBOUNDED,
    },
    {
        .group   = "Counting",
        .name    = "nan-5p",
        .type    = OPT_SIZE,
        .offset  = offsetof(cli_args, pipeline.rate_config.nan_5p),
        .metavar = "N",
        .help    = "write NaN rates and errors for this many bases at the 5' end",
        .label   = "5' Primer Length",
        .minimum = 0,
        .maximum = CLI_UNBOUNDED,
    },
    {
        .group   = "Counting",
        .name    = "nan-3p",
        .type    = OPT_SIZE,
        .offset  = offsetof(cli_args, pipeline.rate_config.nan_3p),
        .metavar = "N",
        .help    = "write NaN rates and errors for this many bases at the 3' end",
        .label   = "3' Primer Length",
        .minimum = 0,
        .maximum = CLI_UNBOUNDED,
    },

    {
        .group       = "Counting",
        .name        = "params",
        .type        = OPT_STRING,
        .offset      = offsetof(cli_args, params_path),
        .metavar     = "FILE",
        .help        = "read the pair HMM's rates from this file",
        .unset_label = "built in",
    },

    {
        .group   = "Performance",
        .name    = "workers",
        .key     = 'j',
        .type    = OPT_SIZE,
        .offset  = offsetof(cli_args, pipeline.workers),
        .metavar = "N",
        .help    = "threads running the processing step",
        .minimum = 1,
        .maximum = CLI_UNBOUNDED,
    },
    {
        .group   = "Performance",
        .name    = "decode-threads",
        .type    = OPT_INT,
        .offset  = offsetof(cli_args, pipeline.decode_threads),
        .metavar = "N",
        .help    = "htslib threads for BGZF decompression",
        .minimum = 0,
        .maximum = CLI_UNBOUNDED,
    },
    {
        .group   = "Performance",
        .name    = "queue-capacity",
        .type    = OPT_SIZE,
        .offset  = offsetof(cli_args, pipeline.queue_capacity),
        .metavar = "N",
        .help    = "reads in transit at once",
        .minimum = 1,
        .maximum = CLI_UNBOUNDED,
        .hidden  = true,
    },
    {
        .group   = "Performance",
        .name    = "batch",
        .type    = OPT_SIZE,
        .offset  = offsetof(cli_args, pipeline.batch),
        .metavar = "N",
        .help    = "reads transferred per queue operation",
        .minimum = 1,
        .maximum = CLI_UNBOUNDED,
        .hidden  = true,
    },
    {
        .group   = "Performance",
        .name    = "live-refs",
        .type    = OPT_SIZE,
        .offset  = offsetof(cli_args, pipeline.live_refs),
        .metavar = "N",
        .help    = "references in flight",
        .minimum = 1,
        .maximum = CLI_UNBOUNDED,
        .hidden  = true,
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
        .name   = "dump-params",
        .type   = OPT_FLAG,
        .help   = "write the rates in the form --params reads and exit",
        .action = CLI_PRINT,
        .print  = params_dump_defaults,
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
        .name         = "alignment",
        .metavar      = "BAM",
        .help         = "coordinate-sorted alignments",
        .offset       = offsetof(cli_args, pipeline.bam_paths),
        .required     = true,
        .variadic     = true,
        .count_offset = offsetof(cli_args, pipeline.n_bams),
    },
};

cli_args cmuts_hmm_defaults(void)
{
    return (cli_args){ .pipeline = pipeline_defaults() };
}

cli_spec cmuts_hmm_spec(const cli_args *defaults)
{
    return (cli_spec){
        .program       = "cmuts hmm",
        .version       = CMUTS_VERSION,
        .summary       = "Fast, multithreaded pair-HMM counting of MaP-seq mutations.",
        .options       = OPTIONS,
        .n_options     = sizeof OPTIONS / sizeof *OPTIONS,
        .positionals   = POSITIONALS,
        .n_positionals = sizeof POSITIONALS / sizeof *POSITIONALS,
        .defaults      = defaults,
        .args_size     = sizeof *defaults,
    };
}
