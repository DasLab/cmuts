/* options.h -- what cmuts align accepts on the command line.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "cli.h"
#include "mapping.h"
#include "version.h"

/* Everything the command line can set.
 *
 * Module configuration is held in a nested struct and not flattened, so that the module
 * keeps owning its own settings and this stays a container. */
typedef struct {
    mapping_config mapping;

    int preset;  /* a preset_id, held as the int the parser writes */
} align_args;

align_args align_defaults(void);

/* The spec borrows defaults, which must outlive it. */
cli_spec align_spec(const align_args *defaults);

/* Completes the module's settings from what the parser filled in elsewhere. */
void align_take_arguments(align_args *args);
