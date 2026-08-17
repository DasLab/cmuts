/* main.c -- driver for generating alignments and the reference they came from.
 *
 * Writes a BAM and a matching FASTA, for tests needing a known shape and benchmarks needing
 * volume.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <stdio.h>

#include <htslib/hts_log.h>

#include "cli.h"
#include "dataset.h"
#include "error.h"
#include "options.h"
#include "subcommands.h"

int gen_main(int argc, char **argv)
{
    gen_args       defaults = gen_defaults();
    cli_spec       spec     = gen_spec(&defaults);
    gen_args       args;
    dataset_config cfg;
    char           error[CM_ERROR_MAX];

    /* Silenced for the reason cmuts hmm silences it: everything htslib reports it also
     * returns. */
    hts_set_log_level(HTS_LOG_OFF);

    switch (cli_parse(&spec, argc, argv, &args)) {
        case CLI_DONE:  return 0;
        case CLI_ERROR: return 2;
        case CLI_OK:    break;
    }

    if (gen_configure(&cfg, &args, error, sizeof error) < 0 ||
        dataset_write(&cfg, error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", spec.program, error);
        return 1;
    }

    return 0;
}
