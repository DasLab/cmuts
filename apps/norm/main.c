/* main.c -- driver for normalization by one norm taken from the rates.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <stdio.h>

#include "cli.h"
#include "error.h"
#include "normalize.h"
#include "options.h"
#include "subcommands.h"

int norm_main(int argc, char **argv)
{
    norm_args defaults = norm_defaults();
    cli_spec  spec     = norm_spec(&defaults);
    norm_args args;
    char      error[CM_ERROR_MAX];

    switch (cli_parse(&spec, argc, argv, &args)) {
        case CLI_DONE:  return 0;
        case CLI_ERROR: return 2;
        case CLI_OK:    break;
    }

    if (norm_take_arguments(&args, error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", spec.program, error);
        return 2;
    }

    if (normalize_run(&args.normalize, spec.program, &CMUTS_NORM_WRITES, error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", spec.program, error);
        return 1;
    }

    return 0;
}
