/* main.c -- driver for read alignment.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <stdio.h>

#include "cli.h"
#include "error.h"
#include "mapping.h"
#include "options.h"
#include "subcommands.h"

int align_main(int argc, char **argv)
{
    align_args defaults = align_defaults();
    cli_spec   spec     = align_spec(&defaults);
    align_args args;
    char       error[CM_ERROR_MAX];

    switch (cli_parse(&spec, argc, argv, &args)) {
        case CLI_DONE:  return 0;
        case CLI_ERROR: return 2;
        case CLI_OK:    break;
    }

    align_take_arguments(&args);

    if (mapping_run(&args.mapping, error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", spec.program, error);
        return 1;
    }

    return 0;
}
