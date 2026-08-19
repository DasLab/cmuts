/* main.c -- driver for scoring against a known structure.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <stdio.h>

#include "cli.h"
#include "error.h"
#include "options.h"
#include "score.h"
#include "subcommands.h"

int score_main(int argc, char **argv)
{
    score_args defaults = score_defaults();
    cli_spec   spec     = score_spec(&defaults);
    score_args args;
    char       error[CM_ERROR_MAX];

    switch (cli_parse(&spec, argc, argv, &args)) {
        case CLI_DONE:  return 0;
        case CLI_ERROR: return 2;
        case CLI_OK:    break;
    }

    if (score_run(&args.score, error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", spec.program, error);
        return 1;
    }

    return 0;
}
