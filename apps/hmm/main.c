/* main.c -- driver for the processing pipeline.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <stdio.h>

#include <htslib/hts_log.h>

#include "cli.h"
#include "error.h"
#include "filter.h"
#include "options.h"
#include "params.h"
#include "pipeline.h"
#include "subcommands.h"

int hmm_main(int argc, char **argv)
{
    cli_args defaults = cmuts_hmm_defaults();
    cli_spec spec     = cmuts_hmm_spec(&defaults);
    cli_args args;
    char     error[CM_ERROR_MAX];

    /* htslib writes its own report of a failure to stderr, which would appear alongside
     * ours. Everything it reports it also returns, so silencing it loses no information. */
    hts_set_log_level(HTS_LOG_OFF);

    switch (cli_parse(&spec, argc, argv, &args)) {
        case CLI_DONE:  return 0;
        case CLI_ERROR: return 2;
        case CLI_OK:    break;
    }

    /* A usage error, reported here: cli_parse checks each option against its own range,
     * and a constraint between two options cannot be written in the table. */
    if (!filter_satisfiable(&args.pipeline.filter_config)) {
        fprintf(stderr, "%s: --max-length is below --min-length\n", spec.program);
        return 2;
    }

    /* The file names the rates it changes, so what is read stands over the defaults the
     * spec filled in. */
    if (args.params_path &&
        params_read(args.params_path, &args.pipeline.tally_config.params,
                    error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", spec.program, error);
        return 2;
    }

    if (pipeline_run(&args.pipeline, spec.program, &CMUTS_HMM_WRITES, error,
                     sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", spec.program, error);
        return 1;
    }

    return 0;
}
