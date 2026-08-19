/* options.h -- what cmuts score accepts on the command line.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "cli.h"
#include "score.h"
#include "version.h"

/* Everything the command line can set.
 *
 * Module configuration is held in a nested struct and not flattened, so that the module
 * keeps owning its own settings and this stays a container. */
typedef struct {
    score_config score;
} score_args;

score_args score_defaults(void);

/* The spec borrows defaults, which must outlive it. */
cli_spec score_spec(const score_args *defaults);
