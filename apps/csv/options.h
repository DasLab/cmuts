/* options.h -- what cmuts csv accepts on the command line.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "cli.h"
#include "csv.h"
#include "version.h"

/* Everything the command line can set.
 *
 * Module configuration is held in a nested struct and not flattened, so that the module
 * keeps owning its own settings and this stays a container. */
typedef struct {
    csv_config csv;
} csv_args;

csv_args csv_defaults(void);

/* The spec borrows defaults, which must outlive it. */
cli_spec csv_spec(const csv_args *defaults);
