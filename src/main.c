/* main.c -- driver for the processing pipeline.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <stdio.h>

#include <htslib/hts_log.h>

#include "cli.h"
#include "error.h"
#include "options.h"
#include "pipeline.h"

int main(int argc, char **argv)
{
    cli_args defaults = cmuts_defaults();
    cli_spec spec     = cmuts_spec(&defaults);
    cli_args args;
    char     error[CM_ERROR_MAX];

    /* htslib writes its own account of a failure to stderr, which would arrive
     * beside ours saying the same thing in another voice. What it reports it
     * also returns, so nothing is lost by keeping it quiet. */
    hts_set_log_level(HTS_LOG_OFF);

    switch (cli_parse(&spec, argc, argv, &args)) {
        case CLI_DONE:  return 0;
        case CLI_ERROR: return 2;
        case CLI_OK:    break;
    }

    if (pipeline_run(&args.pipeline, error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", spec.program, error);
        return 1;
    }

    return 0;
}
