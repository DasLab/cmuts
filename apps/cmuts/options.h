/* options.h -- what cmuts accepts on the command line.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "cli.h"
#include "pipeline.h"
#include "version.h"

/* Everything the command line can set.
 *
 * Module configuration is held by nested structs rather than flattened, so
 * that each module keeps owning its own settings and this stays a container.
 * As processing options arrive they get their own struct here, not more loose
 * fields. */
typedef struct {
    pipeline_config pipeline;
} cli_args;

cli_args cmuts_defaults(void);

/* The spec borrows defaults, which must outlive it. */
cli_spec cmuts_spec(const cli_args *defaults);
