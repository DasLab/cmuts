/* options.c -- cmuts-hmm's command line, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "options.h"

#include <stddef.h>

#include "filter.h"
#include "output.h"
#include "pairs.h"
#include "params.h"
#include "phmm.h"

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

/* The datasets a run leaves behind. Every output holds the first of them; the pairwise
 * squares are written only when they are asked for. */
static const out_written WRITTEN[] = {
    { .id = OUT_COVERAGE,
      .detail = "The number of reads in which this base was present, weighted by PHRED scores." },
    { .id = OUT_REACTIVITY,
      .detail = "The mutation rate at this base, weighted by PHRED scores and in accordance with the HMM parameters." },
    { .id = OUT_ERROR,
      .detail = "Standard error of the reactivity values. Purely the statistical error introduced by finite read depths; does not account for experimental or systemic errors." },
    { .id = OUT_LENGTHS,
      .detail = "The number of reads passing all filters, binned by length." },
    { .id = OUT_READS,
      .detail = "The number of reads passing all filters." },
    { .id = OUT_REJECTED,
      .detail = "The number of reads rejected by at least one filter, or which couldn't be modelled by the HMM" },
    { .id = OUT_UNMAPPED,
      .detail = "The number of reads not aligned to any reference." },
    { .id = OUT_PAIRWISE_CORRELATION,
      .detail = "The Pearson correlation of mutations between this pair of bases.",
      .condition = "--pairwise correlation" },
    { .id = OUT_PAIRWISE_CONDITIONAL,
      .detail = "The probability that the base on the first axis was mutated in a read, given that the base on the second axis was.",
      .condition = "--pairwise conditional" },
    { .id = OUT_PAIRWISE_COVERAGE,
      .detail = "The number of reads in which this pair of bases was present, weighted by PHRED scores.",
      .condition = "--pairwise" },
};


const out_manifest CMUTS_HMM_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };

static void dump_layout(FILE *out)
{
    out_dump_layout(out, "cmuts-hmm", &CMUTS_HMM_WRITES);
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
        .help    = "discard alignments below this mapping quality",
        .minimum = 0,
        .maximum = FILTER_MAPQ_UNAVAILABLE - 1,
    },
    {
        .group       = "Filtering",
        .name        = "min-length",
        .type        = OPT_INT,
        .offset      = offsetof(cli_args, pipeline.filter_config.min_length),
        .metavar     = "N",
        .help        = "discard reads shorter than this",
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
        .help        = "discard reads longer than this",
        .unset_label = "no limit",
        .minimum     = 0,
        .maximum     = CLI_UNBOUNDED,
    },
    {
        .group   = "Filtering",
        .name    = "strand",
        .key     = 's',
        .type    = OPT_SET,
        .offset  = offsetof(cli_args, pipeline.filter_config.strand),
        .metavar = "STRANDS",
        .help    = "keep alignments on these strands",
        .choices = STRAND_CHOICES,
    },

    {
        .group   = "Counting",
        .name    = "band",
        .type    = OPT_INT,
        .offset  = offsetof(cli_args, pipeline.tally_config.band),
        .metavar = "N",
        .help    = "reference positions the marginal may look either side of the CIGAR",
        .minimum = 0,
        .maximum = CLI_UNBOUNDED,
    },
    {
        .group   = "Counting",
        .name    = "pairwise",
        .type    = OPT_SET,
        .offset  = offsetof(cli_args, pipeline.pairwise),
        .metavar = "STATS",
        .help    = "write these statistics of how often two positions are modified together",
        .choices = PAIRWISE_CHOICES,
    },
    {
        .group   = "Counting",
        .name    = "min-depth",
        .type    = OPT_DOUBLE,
        .offset  = offsetof(cli_args, pipeline.rate_config.min_depth),
        .metavar = "D",
        .help    = "evidence a position needs before its rate is written",
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
        .group   = "Counting",
        .name    = "substitution-weight",
        .type    = OPT_DOUBLE,
        .offset  = offsetof(cli_args, pipeline.tally_config.weights.weight[PHMM_SUBSTITUTION]),
        .metavar = "W",
        .help    = "what a substitution counts towards the mutation total",
        .minimum = 0,
        .maximum = 1,
    },
    {
        .group   = "Counting",
        .name    = "deletion-weight",
        .type    = OPT_DOUBLE,
        .offset  = offsetof(cli_args, pipeline.tally_config.weights.weight[PHMM_DELETION]),
        .metavar = "W",
        .help    = "what a deletion counts towards the mutation total",
        .minimum = 0,
        .maximum = 1,
    },
    {
        .group   = "Counting",
        .name    = "insertion-weight",
        .type    = OPT_DOUBLE,
        .offset  = offsetof(cli_args, pipeline.tally_config.weights.weight[PHMM_INSERTION]),
        .metavar = "W",
        .help    = "what an insertion counts towards the mutation total",
        .minimum = 0,
        .maximum = 1,
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
        .program       = "cmuts-hmm",
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
