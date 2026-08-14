/* options.h -- what cmuts-hmm accepts on the command line.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "cli.h"
#include "pipeline.h"
#include "version.h"

/* Everything the command line can set.
 *
 * Module configuration is held in nested structs and not flattened, so that each module
 * keeps owning its own settings and this stays a container. New processing options get
 * their own struct here. */
typedef struct {
    pipeline_config pipeline;
} cli_args;

cli_args cmuts_hmm_defaults(void);

/* The spec borrows defaults, which must outlive it. */
cli_spec cmuts_hmm_spec(const cli_args *defaults);
