/* pipeline.h -- streaming a BAM through parallel processing into per-reference results.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "filter.h"

typedef struct {
    const char   *bam_path;
    const char   *fasta_path;
    const char   *output_path;
    bool          overwrite;
    filter_config filter;
    size_t      workers;         /* threads running the processing step */
    int         decode_threads;  /* htslib threads for BGZF inflation */
    size_t      queue_capacity;  /* reads that may be in transit at once */
    size_t      batch;           /* reads transferred per queue operation */
    size_t      live_refs;       /* references in flight; 0 derives it from the
                                    longest reference and a memory budget */
} pipeline_config;

pipeline_config pipeline_defaults(void);

/* Runs the file to completion, leaving the results in the output file.
 * Returns 0, or -1 with a description in error. */
int pipeline_run(const pipeline_config *cfg, char *error, size_t error_len);
