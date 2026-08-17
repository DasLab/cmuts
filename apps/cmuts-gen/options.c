/* options.c -- cmuts-gen's command line, as one table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "options.h"

#include <stddef.h>

#include "version.h"

static const cli_choice FORMAT_CHOICES[] = {
    { "bam", DATASET_BAM },
    { "sam", DATASET_SAM },
    { NULL,  0           },
};

static const cli_option OPTIONS[] = {
    {
        .group    = "Output",
        .name     = "output",
        .key      = 'o',
        .type     = OPT_STRING,
        .offset   = offsetof(gen_args, output),
        .metavar  = "PREFIX",
        .help     = "write PREFIX.bam and PREFIX.fasta",
        .required = true,
    },
    {
        .group   = "Output",
        .name    = "format",
        .type    = OPT_ENUM,
        .offset  = offsetof(gen_args, format),
        .metavar = "FORMAT",
        .help    = "alignment format to write",
        .choices = FORMAT_CHOICES,
    },

    {
        .group   = "Layout",
        .name    = "references",
        .type    = OPT_SIZE,
        .offset  = offsetof(gen_args, references),
        .metavar = "N",
        .help    = "how many references to write",
        .minimum = 1,
        .maximum = CLI_UNBOUNDED,
    },
    {
        .group   = "Layout",
        .name    = "ref-length",
        .type    = OPT_STRING,
        .offset  = offsetof(gen_args, ref_length),
        .metavar = "DISTRIBUTION",
        .help    = "length of each reference",
    },
    {
        .group   = "Layout",
        .name    = "covered",
        .type    = OPT_DOUBLE,
        .offset  = offsetof(gen_args, covered),
        .metavar = "F",
        .help    = "fraction of references receiving any reads",
        .minimum = 0,
        .maximum = 1,
    },
    {
        .group   = "Layout",
        .name    = "reads-per-ref",
        .type    = OPT_STRING,
        .offset  = offsetof(gen_args, reads_per_ref),
        .metavar = "DISTRIBUTION",
        .help    = "reads on each covered reference",
    },

    {
        .group   = "Reads",
        .name    = "read-length",
        .type    = OPT_STRING,
        .offset  = offsetof(gen_args, read_length),
        .metavar = "DISTRIBUTION",
        .help    = "reference span each read covers",
    },
    {
        .group   = "Reads",
        .name    = "mapq",
        .type    = OPT_STRING,
        .offset  = offsetof(gen_args, mapq),
        .metavar = "DISTRIBUTION",
        .help    = "mapping quality of each read",
    },
    {
        .group   = "Reads",
        .name    = "base-quality",
        .type    = OPT_STRING,
        .offset  = offsetof(gen_args, base_quality),
        .metavar = "DISTRIBUTION",
        .help    = "PHRED score of each base",
    },
    {
        .group   = "Reads",
        .name    = "reverse",
        .type    = OPT_DOUBLE,
        .offset  = offsetof(gen_args, reverse),
        .metavar = "F",
        .help    = "fraction of reads on the reverse strand",
        .minimum = 0,
        .maximum = 1,
    },
    {
        .group   = "Reads",
        .name    = "unmapped",
        .type    = OPT_STRING,
        .offset  = offsetof(gen_args, unmapped),
        .metavar = "DISTRIBUTION",
        .help    = "reads aligning nowhere",
    },

    {
        .group   = "Differences from the reference",
        .name    = "mismatch-rate",
        .type    = OPT_DOUBLE,
        .offset  = offsetof(gen_args, mismatch_rate),
        .metavar = "F",
        .help    = "per aligned base",
        .minimum = 0,
        .maximum = 1,
    },
    {
        .group   = "Differences from the reference",
        .name    = "insertions",
        .type    = OPT_STRING,
        .offset  = offsetof(gen_args, insertions),
        .metavar = "DISTRIBUTION",
        .help    = "insertion events per read",
    },
    {
        .group   = "Differences from the reference",
        .name    = "insertion-length",
        .type    = OPT_STRING,
        .offset  = offsetof(gen_args, insertion_length),
        .metavar = "DISTRIBUTION",
        .help    = "bases per insertion",
    },
    {
        .group   = "Differences from the reference",
        .name    = "deletions",
        .type    = OPT_STRING,
        .offset  = offsetof(gen_args, deletions),
        .metavar = "DISTRIBUTION",
        .help    = "deletion events per read",
    },
    {
        .group   = "Differences from the reference",
        .name    = "deletion-length",
        .type    = OPT_STRING,
        .offset  = offsetof(gen_args, deletion_length),
        .metavar = "DISTRIBUTION",
        .help    = "bases per deletion",
    },
    {
        .group   = "Differences from the reference",
        .name    = "soft-clips",
        .type    = OPT_STRING,
        .offset  = offsetof(gen_args, soft_clips),
        .metavar = "DISTRIBUTION",
        .help    = "clipped ends per read, none through both",
    },
    {
        .group   = "Differences from the reference",
        .name    = "soft-clip-length",
        .type    = OPT_STRING,
        .offset  = offsetof(gen_args, soft_clip_length),
        .metavar = "DISTRIBUTION",
        .help    = "bases per clipped end",
    },

    {
        .group   = "Determinism",
        .name    = "seed",
        .type    = OPT_SIZE,
        .offset  = offsetof(gen_args, seed),
        .metavar = "N",
        .help    = "everything generated follows from this",
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

gen_args gen_defaults(void)
{
    return (gen_args){
        .format           = DATASET_BAM,
        .references       = 100,
        .ref_length       = "400",
        .covered          = 1.0,
        .reads_per_ref    = "20",
        .read_length      = "100:400",
        .mapq             = "0,1,10,30,42,60",
        .base_quality     = "30:40",
        .reverse          = 0.5,
        .unmapped         = "10",
        .mismatch_rate    = 0.01,
        .insertions       = "0:1",
        .insertion_length = "1:5",
        .deletions        = "0:1",
        .deletion_length  = "1:3",
        .soft_clips       = "0:2",
        .soft_clip_length = "5:30",
        .seed             = 1,
    };
}

cli_spec gen_spec(const gen_args *defaults)
{
    return (cli_spec){
        .program   = "cmuts-gen",
        .version   = CMUTS_VERSION,
        .summary   = "generate alignments and the reference they came from",
        .options   = OPTIONS,
        .n_options = sizeof OPTIONS / sizeof *OPTIONS,
        .defaults  = defaults,
        .args_size = sizeof *defaults,
    };
}

/* ------------------------------------------------------------------------ */
/* Reading the specs                                                         */
/* ------------------------------------------------------------------------ */

/* One option carrying a distribution, paired with the field it is parsed into. The two offsets
 * are declared together so that an option added here is read by the same loop as the rest,
 * and cannot be given a table row and then left unparsed. */
typedef struct {
    const char *name;    /* the option it came from, for the message */
    size_t      text;    /* the written form, within gen_args */
    size_t      parsed;  /* what it becomes, within dataset_config */
} distribution_option;

static const distribution_option DISTRIBUTION_OPTIONS[] = {
    { "ref-length",       offsetof(gen_args, ref_length),
                          offsetof(dataset_config, ref_length)          },
    { "reads-per-ref",    offsetof(gen_args, reads_per_ref),
                          offsetof(dataset_config, reads_per_ref)       },
    { "unmapped",         offsetof(gen_args, unmapped),
                          offsetof(dataset_config, unmapped)            },
    { "read-length",      offsetof(gen_args, read_length),
                          offsetof(dataset_config, model.length)        },
    { "mapq",             offsetof(gen_args, mapq),
                          offsetof(dataset_config, model.mapq)          },
    { "base-quality",     offsetof(gen_args, base_quality),
                          offsetof(dataset_config, model.base_quality)  },
    { "insertions",       offsetof(gen_args, insertions),
                          offsetof(dataset_config, model.insertions)    },
    { "insertion-length", offsetof(gen_args, insertion_length),
                          offsetof(dataset_config, model.insertion_length) },
    { "deletions",        offsetof(gen_args, deletions),
                          offsetof(dataset_config, model.deletions)     },
    { "deletion-length",  offsetof(gen_args, deletion_length),
                          offsetof(dataset_config, model.deletion_length)  },
    { "soft-clips",       offsetof(gen_args, soft_clips),
                          offsetof(dataset_config, model.soft_clips)    },
    { "soft-clip-length", offsetof(gen_args, soft_clip_length),
                          offsetof(dataset_config, model.soft_clip_length) },
};

/* Parses one option's text into its distribution. The option's name is written into error
 * first and the reason appended after it, so no second buffer is needed. */
static int parse_distribution(distribution *out,
                              const distribution_option *option,
                              const char *text,
                              char *error, size_t error_len)
{
    int prefix = snprintf(error, error_len, "--%s: ", option->name);

    if (prefix < 0 || (size_t)prefix >= error_len) {
        return -1;
    }

    if (distribution_parse(out, text, error + prefix,
                           error_len - (size_t)prefix) != 0) {
        return -1;
    }

    error[0] = '\0';
    return 0;
}

static const char *written_form(const gen_args *args,
                                const distribution_option *option)
{
    return *(const char *const *)((const char *)args + option->text);
}

static distribution *destination(dataset_config *cfg,
                                 const distribution_option *option)
{
    return (distribution *)((char *)cfg + option->parsed);
}

int gen_configure(dataset_config *cfg, const gen_args *args,
                  char *error, size_t error_len)
{
    *cfg = (dataset_config){
        .prefix                 = args->output,
        .format                 = (dataset_format)args->format,
        .references             = args->references,
        .covered                = args->covered,
        .seed                   = args->seed,
        .model.mismatch_rate    = args->mismatch_rate,
        .model.reverse_fraction = args->reverse,
    };

    size_t n = sizeof DISTRIBUTION_OPTIONS / sizeof *DISTRIBUTION_OPTIONS;

    for (size_t i = 0; i < n; i++) {
        const distribution_option *option = &DISTRIBUTION_OPTIONS[i];

        if (parse_distribution(destination(cfg, option), option,
                       written_form(args, option), error, error_len) < 0) {
            return -1;
        }
    }

    return 0;
}
