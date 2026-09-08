/* main.c -- driver for writing an output as comma separated values.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <stdio.h>

#include "cli.h"
#include "csv.h"
#include "error.h"
#include "manifest.h"
#include "options.h"
#include "subcommands.h"

int csv_main(int argc, char **argv)
{
    csv_args defaults = csv_defaults();
    cli_spec spec     = csv_spec(&defaults);
    csv_args args;
    char     error[CM_ERROR_MAX];

    switch (cli_parse(&spec, argc, argv, &args)) {
        case CLI_DONE:  return 0;
        case CLI_ERROR: return 2;
        case CLI_OK:    break;
    }

    if (csv_run(&args.csv, CMUTS_CSV_READS, CMUTS_CSV_N_READS, stdout,
                error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", spec.program, error);
        return 1;
    }

    return 0;
}
