/* pipeline.h -- streaming a BAM through parallel processing into per-reference results.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include "filter.h"

typedef struct {
    const char   *bam_path;
    const char   *fasta_path;
    const char   *output_path;
    filter_config filter;
    size_t      workers;         /* threads running the processing step */
    int         decode_threads;  /* htslib threads for BGZF inflation */
    size_t      queue_capacity;  /* reads that may be in transit at once */
    size_t      batch;           /* reads transferred per queue operation */
    size_t      live_refs;       /* references in flight; 0 derives it from the
                                    longest reference and a memory budget */
} pipeline_config;

pipeline_config pipeline_defaults(void);

typedef struct {
    size_t reads_total;
    size_t reads_unmapped;
    size_t reads_filtered;
    size_t reads_processed;
    size_t refs_completed;

    /* Sums over every accumulator written. Each is a total of integer counts
     * held as double, so it is exact and independent of the order in which
     * threads happened to merge: two runs at different worker counts must
     * agree bit for bit. */
    double coverage_total;
    double mutations_total;
    double reads_recorded;
} pipeline_stats;

/* Runs the file to completion. Returns 0, or -1 with a description in error. */
int pipeline_run(const pipeline_config *cfg, pipeline_stats *stats,
                 char *error, size_t error_len);
