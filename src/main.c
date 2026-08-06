/* main.c -- driver for the processing pipeline.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <stdio.h>

#include "cli.h"
#include "pipeline.h"

#define ERROR_MAX 512

int main(int argc, char **argv)
{
    cli_args args;
    char     error[ERROR_MAX];

    switch (cli_parse(argc, argv, &args)) {
        case CLI_DONE:  return 0;
        case CLI_ERROR: return 2;
        case CLI_OK:    break;
    }

    if (pipeline_run(&args.pipeline, error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", argv[0], error);
        return 1;
    }

    return 0;
}
