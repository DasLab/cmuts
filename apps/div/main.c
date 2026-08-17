/* main.c -- driver for normalization against a denatured control.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <stdio.h>

#include "cli.h"
#include "divide.h"
#include "error.h"
#include "options.h"
#include "subcommands.h"

int div_main(int argc, char **argv)
{
    div_args defaults = div_defaults();
    cli_spec spec     = div_spec(&defaults);
    div_args args;
    char     error[CM_ERROR_MAX];

    switch (cli_parse(&spec, argc, argv, &args)) {
        case CLI_DONE:  return 0;
        case CLI_ERROR: return 2;
        case CLI_OK:    break;
    }

    if (divide_run(&args.divide, spec.program, &CMUTS_DIV_WRITES, error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", spec.program, error);
        return 1;
    }

    return 0;
}
