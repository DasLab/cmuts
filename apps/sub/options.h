/* options.h -- what cmuts sub accepts on the command line.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "output.h"

#include "cli.h"
#include "subtract.h"
#include "version.h"

/* Everything the command line can set.
 *
 * Module configuration is held in a nested struct and not flattened, so that the module
 * keeps owning its own settings and this stays a container. */
typedef struct {
    subtract_config subtract;
} sub_args;

sub_args sub_defaults(void);

/* The spec borrows defaults, which must outlive it. */
cli_spec sub_spec(const sub_args *defaults);

/* The datasets a run of this program writes. */
extern const out_manifest CMUTS_SUB_WRITES;
