/* pipeline.h -- streaming a BAM through parallel processing into per-reference results.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "filter.h"
#include "rates.h"
#include "refseq.h"
#include "tally.h"

typedef struct {
    /* Read as one file: every reference is counted across all of them, into the single row
     * it has in the output. */
    const char *const *bam_paths;
    size_t             n_bams;
    const char        *fasta_path;
    const char        *output_path;
    bool               overwrite;
    int                verify;  /* refseq_check bits the FASTA is held to */
    rate_config        rate_config;
    filter_config      filter_config;
    tally_config       tally_config;
    size_t      workers;         /* threads running the processing step */
    int         decode_threads;  /* htslib threads for BGZF inflation */
    size_t      queue_capacity;  /* reads that may be in transit at once */
    size_t      batch;           /* reads transferred per queue operation */
    size_t      live_refs;       /* references in flight at once */
    bool        pairwise;        /* also count how often two positions are modified
                                    together, which costs the square of the longest
                                    reference in every accumulator */
} pipeline_config;

pipeline_config pipeline_defaults(void);

/* Runs the input to completion, leaving the results in the output file. program is
 * recorded in the output as what produced it, and comes from the command line the caller
 * declared. Returns 0, or -1 with a description in error. */
int pipeline_run(const pipeline_config *cfg, const char *program, char *error,
                 size_t error_len);
