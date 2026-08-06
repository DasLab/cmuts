/* options.c -- cmuts's command line, as one table.
 *
 * Rows use designated initializers so that a field added to cli_option
 * defaults quietly rather than having to be spelled out in every one.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "options.h"

#include <stddef.h>
#include <stdint.h>

#include "filter.h"

static const cli_choice STRAND_CHOICES[] = {
    { "both",    FILTER_STRAND_BOTH    },
    { "forward", FILTER_STRAND_FORWARD },
    { "reverse", FILTER_STRAND_REVERSE },
    { NULL,      0                     },
};

static const cli_option OPTIONS[] = {
    { .group = "Input and output", .name = "fasta", .key = 'f', .type = OPT_STRING,
      .offset = offsetof(cli_args, pipeline.fasta_path), .metavar = "FASTA",
      .help = "reference sequences",
      .detail = "Reference sequences, in the same order as the alignment file's "
                "header. Read as a single forward pass and never indexed, so the "
                "file may be of any size.",
      .required = true },
    { .group = "Input and output", .name = "output", .key = 'o', .type = OPT_STRING,
      .offset = offsetof(cli_args, pipeline.output_path), .metavar = "HDF5",
      .help = "write results to this file",
      .detail = "Results are written as one row per reference, with references "
                "that received no reads left as NaN.",
      .required = true },

    { .group = "Filtering", .name = "min-mapq", .key = 'q', .type = OPT_INT,
      .offset = offsetof(cli_args, pipeline.filter.min_mapq), .metavar = "N",
      .help = "discard alignments below this mapping quality",
      .detail = "Compared numerically, as samtools does. MAPQ 255 means "
                "\"unavailable\" rather than \"perfect\", but is treated as passing "
                "any threshold, since an aligner that emits it throughout would "
                "otherwise have all of its output discarded. Unmapped reads are "
                "excluded regardless and counted separately.",
      .minimum = 0, .maximum = 255 },
    { .group = "Filtering", .name = "min-length", .type = OPT_INT,
      .offset = offsetof(cli_args, pipeline.filter.min_length), .metavar = "N",
      .help = "discard reads shorter than this",
      .detail = "Measured on the stored sequence, so a hard-clipped read counts "
                "only the bases the aligner kept. Left unset, no lower bound is "
                "applied.",
      .unset_label = "no limit", .minimum = 0, .maximum = INT32_MAX },
    { .group = "Filtering", .name = "max-length", .type = OPT_INT,
      .offset = offsetof(cli_args, pipeline.filter.max_length), .metavar = "N",
      .help = "discard reads longer than this",
      .detail = "Measured on the stored sequence, as with --min-length. Left "
                "unset, no upper bound is applied.",
      .unset_label = "no limit", .minimum = 0, .maximum = INT32_MAX },
    { .group = "Filtering", .name = "strand", .key = 's', .type = OPT_ENUM,
      .offset = offsetof(cli_args, pipeline.filter.strand), .metavar = "STRAND",
      .help = "keep alignments on this strand",
      .detail = "Tests the alignment's own reverse bit, which for single-end "
                "reads is the strand the read came from. It says nothing about "
                "the fragment: with paired data, which strand a fragment belongs "
                "to depends on the library protocol and on which mate is in hand.",
      .choices = STRAND_CHOICES },

    { .group = "Performance", .name = "workers", .key = 'j', .type = OPT_SIZE,
      .offset = offsetof(cli_args, pipeline.workers), .metavar = "N",
      .help = "threads running the processing step",
      .detail = "Reads are taken from a shared pool, so a worker is free to cross "
                "reference boundaries and no thread idles waiting for a reference "
                "of its own.",
      .minimum = 1, .maximum = 1024 },
    { .group = "Performance", .name = "decode-threads", .type = OPT_INT,
      .offset = offsetof(cli_args, pipeline.decode_threads), .metavar = "N",
      .help = "htslib threads for BGZF decompression",
      .detail = "Parallelizes inflation only; reading and record parsing stay "
                "sequential. Worth raising when the loader is the bottleneck, and "
                "pointless on small files.",
      .minimum = 0, .maximum = 64 },
    { .group = "Performance", .name = "queue-capacity", .type = OPT_SIZE,
      .offset = offsetof(cli_args, pipeline.queue_capacity), .metavar = "N",
      .help = "reads in transit at once",
      .detail = "Bounds how far the loader may run ahead of the workers, and with "
                "it how much memory reads in flight occupy.",
      .minimum = 1, .maximum = 1 << 20 },
    { .group = "Performance", .name = "batch", .type = OPT_SIZE,
      .offset = offsetof(cli_args, pipeline.batch), .metavar = "N",
      .help = "reads transferred per queue operation",
      .detail = "Larger batches amortize locking across more reads, at the cost "
                "of holding that many reads behind a slow one.",
      .minimum = 1, .maximum = 1 << 16 },
    { .group = "Performance", .name = "live-refs", .type = OPT_SIZE,
      .offset = offsetof(cli_args, pipeline.live_refs), .metavar = "N",
      .help = "references in flight",
      .detail = "How far the loader may run ahead of a worker that stalls on one "
                "read. Left unset, a count is derived from the longest reference "
                "and a memory budget, which keeps many short references generous "
                "without letting a few very long ones exhaust memory.",
      .unset_label = "derived from memory", .minimum = 0, .maximum = 1 << 16 },

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

static const cli_positional POSITIONALS[] = {
    { .name = "alignment", .metavar = "BAM", .help = "coordinate-sorted alignments",
      .detail = "Read once, sequentially. Must be coordinate sorted, so that a "
                "reference is finished as soon as the reader moves past it.",
      .offset = offsetof(cli_args, pipeline.bam_path), .required = true },
};

cli_args cmuts_defaults(void)
{
    return (cli_args){ .pipeline = pipeline_defaults() };
}

cli_spec cmuts_spec(const cli_args *defaults)
{
    return (cli_spec){
        .program       = "cmuts",
        .version       = CMUTS_VERSION,
        .summary       = "per-reference accumulation over a BAM file",
        .options       = OPTIONS,
        .n_options     = sizeof OPTIONS / sizeof *OPTIONS,
        .positionals   = POSITIONALS,
        .n_positionals = sizeof POSITIONALS / sizeof *POSITIONALS,
        .defaults      = defaults,
        .args_size     = sizeof *defaults,
    };
}
