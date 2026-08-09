/* options.h -- what cmuts-gen accepts on the command line.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include "cli.h"
#include "dataset.h"

/* Everything the command line can set.
 *
 * Every parameter that varies arrives as a spec -- a constant, a range, or a
 * list -- so the same option pins a value down for a test and spreads it for a
 * benchmark. They are text here and are turned into the dataset's model once,
 * after the command line has been read.
 */
typedef struct {
    const char *output;
    int         format;

    size_t      references;
    const char *ref_length;
    double      covered;
    const char *reads_per_ref;

    const char *read_length;
    double      mismatch_rate;
    const char *insertions;
    const char *insertion_length;
    const char *deletions;
    const char *deletion_length;
    const char *soft_clips;
    const char *soft_clip_length;
    const char *mapq;
    const char *base_quality;
    double      reverse;
    const char *unmapped;

    size_t      seed;
} gen_args;

gen_args gen_defaults(void);

/* The spec borrows defaults, which must outlive it. */
cli_spec gen_spec(const gen_args *defaults);

/* Reads the written specs into the dataset they describe. Returns 0, or -1 with
 * a description in error. */
int gen_configure(dataset_config *cfg, const gen_args *args,
                  char *error, size_t error_len);
