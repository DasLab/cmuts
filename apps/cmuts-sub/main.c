/* main.c -- driver for background subtraction.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <stdio.h>

#include "cli.h"
#include "error.h"
#include "options.h"
#include "subtract.h"

int main(int argc, char **argv)
{
    sub_args defaults = sub_defaults();
    cli_spec spec     = sub_spec(&defaults);
    sub_args args;
    char     error[CM_ERROR_MAX];

    switch (cli_parse(&spec, argc, argv, &args)) {
        case CLI_DONE:  return 0;
        case CLI_ERROR: return 2;
        case CLI_OK:    break;
    }

    if (subtract_run(&args.subtract, error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", spec.program, error);
        return 1;
    }

    return 0;
}
