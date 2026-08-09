/* pipeline.h -- streaming a BAM through parallel processing into per-reference results.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "filter.h"
#include "tally.h"

typedef struct {
    /* Read as though they were one file: every reference is counted across all
     * of them, into the single row it has in the output. */
    const char *const *bam_paths;
    size_t             n_bams;
    const char        *fasta_path;
    const char        *output_path;
    bool               overwrite;
    filter_config      filter_config;
    tally_config       tally_config;
    size_t      workers;         /* threads running the processing step */
    int         decode_threads;  /* htslib threads for BGZF inflation */
    size_t      queue_capacity;  /* reads that may be in transit at once */
    size_t      batch;           /* reads transferred per queue operation */
    size_t      live_refs;       /* references in flight at once */
} pipeline_config;

pipeline_config pipeline_defaults(void);

/* Runs the file to completion, leaving the results in the output file.
 * Returns 0, or -1 with a description in error. */
int pipeline_run(const pipeline_config *cfg, char *error, size_t error_len);
