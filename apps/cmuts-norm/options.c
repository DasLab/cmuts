/* options.c -- cmuts-norm's command line, as one table.
 * Rows use designated initializers, so that a field added to cli_option takes its
 * default without being spelled out in every row.
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 * */

#include "options.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>

/* The coverage a position needs before its rate joins the ubr pool. */
#define DEFAULT_MIN_COVERAGE 500

static const cli_choice SCHEME_CHOICES[] = {
    { "ubr",     NORM_UBR     },
    { "outlier", NORM_OUTLIER },
    { NULL,      0            },
};

/* The datasets a run leaves behind: what it was given, divided by one scale, and that
 * scale alongside it. */
static const out_written WRITTEN[] = {
    { .id = OUT_COVERAGE,
      .detail = "The number of reads in which this base was present, weighted by PHRED scores." },
    { .id = OUT_REACTIVITY,
      .detail = "The mutation rate divided by the scale this file records, so a rate reads against the scale rather than as a raw frequency." },
    { .id = OUT_ERROR,
      .detail = "Standard error of the reactivity values. Purely the statistical error introduced by finite read depths; does not account for experimental or systemic errors." },
    { .id = OUT_LENGTHS,
      .detail = "The number of reads passing all filters, binned by length." },
    { .id = OUT_READS,
      .detail = "The number of reads passing all filters." },
    { .id = OUT_REJECTED,
      .detail = "The number of reads that contributed nothing: rejected by at least one filter, or carrying something the pair HMM's rates give no alignment of." },
    { .id = OUT_UNMAPPED,
      .detail = "The number of reads not aligned to any reference." },
    { .id = OUT_NORM,
      .detail = "The scale every rate in this file was divided by." },
};


const out_manifest CMUTS_NORM_WRITES = { WRITTEN, sizeof WRITTEN / sizeof *WRITTEN };

static void dump_layout(FILE *out)
{
    out_dump_layout(out, "cmuts-norm", &CMUTS_NORM_WRITES);
}

static const cli_option OPTIONS[] = {
    {
        .group        = "Input and output",
        .name         = "output",
        .key          = 'o',
        .type         = OPT_STRING,
        .offset       = offsetof(norm_args, output),
        .metavar      = "HDF5",
        .help         = "write one input's results to this file; repeat once per input",
        .required     = true,
        .repeatable   = true,
        .count_offset = offsetof(norm_args, n_outputs),
        .capacity     = NORM_MAX_FILES,
    },
    {
        .group  = "Input and output",
        .name   = "overwrite",
        .type   = OPT_FLAG,
        .offset = offsetof(norm_args, normalize.overwrite),
        .help   = "replace the output files if they already exist",
    },

    {
        .group   = "Normalization",
        .name    = "norm",
        .type    = OPT_ENUM,
        .offset  = offsetof(norm_args, scheme),
        .metavar = "SCHEME",
        .help    = "how the scale is taken from the rates",
        .choices = SCHEME_CHOICES,
    },
    {
        .group   = "Normalization",
        .name    = "min-coverage",
        .type    = OPT_DOUBLE,
        .offset  = offsetof(norm_args, normalize.min_coverage),
        .metavar = "N",
        .help    = "coverage a position needs before its rate sets the scale (ubr only)",
        .minimum = 0,
        .maximum = CLI_UNBOUNDED,
    },

    {
        .group       = "Clipping",
        .name        = "clip-below",
        .type        = OPT_DOUBLE,
        .offset      = offsetof(norm_args, normalize.clip_below),
        .metavar     = "N",
        .help        = "raise a normalized reactivity up to this value",
        .unset_label = "none",
        .minimum     = -CLI_UNBOUNDED,
        .maximum     = CLI_UNBOUNDED,
    },
    {
        .group       = "Clipping",
        .name        = "clip-above",
        .type        = OPT_DOUBLE,
        .offset      = offsetof(norm_args, normalize.clip_above),
        .metavar     = "N",
        .help        = "lower a normalized reactivity down to this value",
        .unset_label = "none",
        .minimum     = -CLI_UNBOUNDED,
        .maximum     = CLI_UNBOUNDED,
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
        .help   = "describe the output format as JSON and exit",
        .hidden = true,
        .action = CLI_PRINT,
        .print  = dump_layout,
    },
};

static const cli_positional POSITIONALS[] = {
    {
        .name         = "rates",
        .metavar      = "RATES",
        .help         = "the reactivities to normalize",
        .offset       = offsetof(norm_args, normalize.inputs),
        .required     = true,
        .variadic     = true,
        .count_offset = offsetof(norm_args, normalize.n_files),
    },
};

norm_args norm_defaults(void)
{
    return (norm_args){
        .normalize = {
            .min_coverage = DEFAULT_MIN_COVERAGE,
            .clip_below   = (double)NAN,
            .clip_above   = (double)NAN,
        },
        .scheme = NORM_UBR,
    };
}

cli_spec norm_spec(const norm_args *defaults)
{
    return (cli_spec){
        .program       = "cmuts-norm",
        .version       = CMUTS_VERSION,
        .summary       = "divide cmuts outputs by one scale taken from their rates.",
        .options       = OPTIONS,
        .n_options     = sizeof OPTIONS / sizeof *OPTIONS,
        .positionals   = POSITIONALS,
        .n_positionals = sizeof POSITIONALS / sizeof *POSITIONALS,
        .defaults      = defaults,
        .args_size     = sizeof *defaults,
    };
}

/* The count is checked before the pairing, so that a run of too many files is told which
 * limit it passed rather than that its two counts disagree. */
static int check_counts(const norm_args *args, char *error, size_t error_len)
{
    if (args->normalize.n_files > NORM_MAX_FILES) {
        snprintf(error, error_len, "at most %d inputs may be given", NORM_MAX_FILES);
        return -1;
    }

    if (args->n_outputs != args->normalize.n_files) {
        snprintf(error, error_len,
                 "%zu input%s given and %zu output%s; --output takes one per input",
                 args->normalize.n_files, args->normalize.n_files == 1 ? "" : "s",
                 args->n_outputs, args->n_outputs == 1 ? "" : "s");
        return -1;
    }

    return 0;
}

/* The parser fills the outputs and the scheme outside the module's own settings, having
 * no place to write an array or an enum into. */
static void take_parsed(norm_args *args)
{
    args->normalize.outputs = args->output;
    args->normalize.scheme  = (norm_scheme)args->scheme;
}

int norm_take_arguments(norm_args *args, char *error, size_t error_len)
{
    if (check_counts(args, error, error_len) < 0) {
        return -1;
    }

    take_parsed(args);
    return 0;
}
