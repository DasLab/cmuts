/* options.h -- what cmuts-div accepts on the command line.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "cli.h"
#include "divide.h"
#include "version.h"

/* Everything the command line can set.
 *
 * Module configuration is held in a nested struct and not flattened, so that the module
 * keeps owning its own settings and this stays a container. */
typedef struct {
    divide_config divide;
} div_args;

div_args div_defaults(void);

/* The spec borrows defaults, which must outlive it. */
cli_spec div_spec(const div_args *defaults);
