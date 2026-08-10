/* options.c -- cmuts's command line, as one table.
 *
 * Rows use designated initializers so that a field added to cli_option
 * defaults quietly rather than having to be spelled out in every one.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "options.h"

#include <stddef.h>

#include "filter.h"
#include "phmm.h"

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

    { .group = "Input and output", .name = "overwrite", .type = OPT_FLAG,
      .offset = offsetof(cli_args, pipeline.overwrite),
      .help = "replace the output file if it already exists",
      .detail = "Without this a run stops rather than destroying a result that "
                "may have cost a great deal more to produce than the one being "
                "started." },

    { .group = "Filtering", .name = "min-mapq", .key = 'q', .type = OPT_INT,
      .offset = offsetof(cli_args, pipeline.filter_config.min_mapq), .metavar = "N",
      .help = "discard alignments below this mapping quality",
      .detail = "Compared numerically, as samtools does. MAPQ 255 means "
                "\"unavailable\" rather than \"perfect\", but is treated as passing "
                "any threshold, since an aligner that emits it throughout would "
                "otherwise have all of its output discarded. Unmapped reads are "
                "excluded regardless and counted separately.",
      .minimum = 0, .maximum = 255 },
    { .group = "Filtering", .name = "min-length", .type = OPT_INT,
      .offset = offsetof(cli_args, pipeline.filter_config.min_length), .metavar = "N",
      .help = "discard reads shorter than this",
      .detail = "Measured on the stored sequence, which is the length of the "
                "molecule that was sequenced rather than the span it aligns to. "
                "Inserted and soft-clipped bases count; hard-clipped ones cannot, "
                "being absent from the record. Left unset, no lower bound is "
                "applied.",
      .unset_label = "no limit", .minimum = 0, .maximum = CLI_UNBOUNDED },
    { .group = "Filtering", .name = "max-length", .type = OPT_INT,
      .offset = offsetof(cli_args, pipeline.filter_config.max_length), .metavar = "N",
      .help = "discard reads longer than this",
      .detail = "Measured on the stored sequence, as with --min-length. Since "
                "that counts insertions, an upper bound also removes reads far "
                "longer than the reference they align to: a read carrying a large "
                "insertion is long even where its aligned span is ordinary. Left "
                "unset, no upper bound is applied.",
      .unset_label = "no limit", .minimum = 0, .maximum = CLI_UNBOUNDED },
    { .group = "Filtering", .name = "strand", .key = 's', .type = OPT_ENUM,
      .offset = offsetof(cli_args, pipeline.filter_config.strand), .metavar = "STRAND",
      .help = "keep alignments on this strand",
      .detail = "Tests the alignment's own reverse bit, which for single-end "
                "reads is the strand the read came from. It says nothing about "
                "the fragment: with paired data, which strand a fragment belongs "
                "to depends on the library protocol and on which mate is in hand.",
      .choices = STRAND_CHOICES },

    { .group = "Counting", .name = "band", .type = OPT_INT,
      .offset = offsetof(cli_args, pipeline.tally_config.band), .metavar = "N",
      .help = "reference positions the marginal may look either side of the CIGAR",
      .detail = "A read carrying an indel is counted over every alignment a band "
                "this wide around its CIGAR admits, rather than over the one path "
                "the aligner chose to report. The band follows that path rather "
                "than a diagonal, so what it bounds is the departure from the "
                "alignment already found. Cost is linear in it, and 0 pins the "
                "marginal to the CIGAR, leaving every read with an indel counted "
                "as written. Nothing bounds it above, and a band wide enough to "
                "exhaust memory ends the run.",
      .minimum = 0, .maximum = CLI_UNBOUNDED },
    { .group = "Counting", .name = "substitution-weight", .type = OPT_DOUBLE,
      .offset = offsetof(cli_args,
                         pipeline.tally_config.weights.weight[PHMM_SUBSTITUTION]),
      .metavar = "W",
      .help = "what a substitution counts towards the mutation total",
      .detail = "Scales what a substitution contributes to the mutations written. "
                "Absolute rather than relative to the other two: raising every "
                "weight raises the total. Coverage, the positions spanned, and "
                "the alignment are unaffected.",
      .minimum = 0, .maximum = CLI_UNBOUNDED },
    { .group = "Counting", .name = "deletion-weight", .type = OPT_DOUBLE,
      .offset = offsetof(cli_args,
                         pipeline.tally_config.weights.weight[PHMM_DELETION]),
      .metavar = "W",
      .help = "what a deletion counts towards the mutation total",
      .detail = "Scales what a deletion contributes to the mutations written, once "
                "per deleted run rather than once per base skipped. Otherwise as "
                "--substitution-weight.",
      .minimum = 0, .maximum = CLI_UNBOUNDED },
    { .group = "Counting", .name = "insertion-weight", .type = OPT_DOUBLE,
      .offset = offsetof(cli_args,
                         pipeline.tally_config.weights.weight[PHMM_INSERTION]),
      .metavar = "W",
      .help = "what an insertion counts towards the mutation total",
      .detail = "Scales what an insertion contributes to the mutations written, "
                "once per inserted run. Otherwise as --substitution-weight. The "
                "default of 0 leaves insertions out of the total.",
      .minimum = 0, .maximum = CLI_UNBOUNDED },

    { .group = "Performance", .name = "workers", .key = 'j', .type = OPT_SIZE,
      .offset = offsetof(cli_args, pipeline.workers), .metavar = "N",
      .help = "threads running the processing step",
      .detail = "Reads are taken from a shared pool, so a worker is free to cross "
                "reference boundaries and no thread idles waiting for a reference "
                "of its own.",
      .minimum = 1, .maximum = CLI_UNBOUNDED },
    { .group = "Performance", .name = "decode-threads", .type = OPT_INT,
      .offset = offsetof(cli_args, pipeline.decode_threads), .metavar = "N",
      .help = "htslib threads for BGZF decompression",
      .detail = "Throughput is flat above the default, so this is a cap rather "
                "than a dial: lower it when the cores are not there to spare. A "
                "total however many files are given, which share the threads "
                "between them.",
      .minimum = 0, .maximum = CLI_UNBOUNDED },
    { .group = "Performance", .name = "queue-capacity", .type = OPT_SIZE,
      .offset = offsetof(cli_args, pipeline.queue_capacity), .metavar = "N",
      .help = "reads in transit at once",
      .detail = "Bounds how far the loader may run ahead of the workers, and with "
                "it how much memory reads in flight occupy.",
      .minimum = 1, .maximum = CLI_UNBOUNDED },
    { .group = "Performance", .name = "batch", .type = OPT_SIZE,
      .offset = offsetof(cli_args, pipeline.batch), .metavar = "N",
      .help = "reads transferred per queue operation",
      .detail = "Larger batches amortize locking across more reads, at the cost "
                "of holding that many reads behind a slow one.",
      .minimum = 1, .maximum = CLI_UNBOUNDED },
    { .group = "Performance", .name = "live-refs", .type = OPT_SIZE,
      .offset = offsetof(cli_args, pipeline.live_refs), .metavar = "N",
      .help = "references in flight",
      .detail = "How far the loader may run ahead of a worker that stalls on one "
                "read. Each reference in flight costs memory in proportion to the "
                "longest reference in the file, so lower it where those are very "
                "long.",
      .minimum = 1, .maximum = CLI_UNBOUNDED },

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
                "reference is finished as soon as the reader moves past it. "
                "Several files are read as though they were one: each reference "
                "is counted across all of them into the single row it has in the "
                "output, which requires that they declare the same references in "
                "the same order.",
      .offset = offsetof(cli_args, pipeline.bam_paths), .required = true,
      .variadic = true, .count_offset = offsetof(cli_args, pipeline.n_bams) },
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
        .summary       = "Fast, multithreaded pair-HMM counting of MaP-seq mutations.",
        .options       = OPTIONS,
        .n_options     = sizeof OPTIONS / sizeof *OPTIONS,
        .positionals   = POSITIONALS,
        .n_positionals = sizeof POSITIONALS / sizeof *POSITIONALS,
        .defaults      = defaults,
        .args_size     = sizeof *defaults,
    };
}
