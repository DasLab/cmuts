/* options.h -- what cmuts-norm accepts on the command line.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "output.h"

#include "cli.h"
#include "normalize.h"
#include "version.h"

/* How many inputs one run takes, which is also how many outputs it writes. A run past this
 * is refused rather than truncated. */
#define NORM_MAX_FILES 256

/* Everything the command line can set.
 *
 * The paths arrive as two arrays the parser fills, so they are held here and pointed at
 * from the module's own settings once both counts are known.
 *
 * Module configuration is held in a nested struct and not flattened, so that the module
 * keeps owning its own settings and this stays a container. */
typedef struct {
    normalize_config normalize;

    const char *output[NORM_MAX_FILES];
    size_t      n_outputs;
    int         scheme;  /* a norm_scheme, held as the int the parser writes */
} norm_args;

norm_args norm_defaults(void);

/* The spec borrows defaults, which must outlive it. */
cli_spec norm_spec(const norm_args *defaults);

/* Completes the module's settings from what the parser filled in elsewhere, refusing a run
 * whose inputs and outputs do not pair up. Returns 0, or -1 with a description in
 * error. */
int norm_take_arguments(norm_args *args, char *error, size_t error_len);

/* The datasets a run of this program writes. */
extern const out_manifest CMUTS_NORM_WRITES;
