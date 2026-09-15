/* mapping.h -- aligning sequencing reads against a reference library.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "preset.h"

/* The largest number of files of reads that one run accepts, which is one file per mate. */
#define MAPPING_MAX_READS 2

/* reads_paths[1] is NULL where the reads arrive as one file. */
typedef struct {
    const char *reads_paths[MAPPING_MAX_READS];
    const char *fasta_path;
    const char *output_path;
    preset_id   preset;
    int         threads;
    size_t      memory;       /* mebibytes the sorter may hold per thread */
    bool        overwrite;
} mapping_config;

/* Aligns the reads against the reference library and writes the sorted alignments to the
 * output file. minimap2 and samtools must be on PATH, and vsearch as well when the reads
 * arrive as a pair of files.
 *
 * One process may run one alignment at a time, because the file that holds the alignments
 * until the run finishes is removed by a signal handler. Returns 0, or -1 with a
 * description in error. */
int mapping_run(const mapping_config *cfg, char *error, size_t error_len);
