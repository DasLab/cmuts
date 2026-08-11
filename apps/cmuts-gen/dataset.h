/* dataset.h -- alignments written beside the reference they came from.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include "sample.h"
#include "simulate.h"

/* SAM is readable in a diff, which is worth having for a small fixture. BAM is for anything
 * large. */
typedef enum {
    DATASET_BAM,
    DATASET_SAM,
} dataset_format;

/* Everything a dataset is made of.
 *
 * How a read differs from the reference is the simulator's to describe, and is held as its
 * model rather than flattened in here. What is left is the shape of the run: how many
 * references there are, how long they are, and how many reads each receives. */
typedef struct {
    const char    *prefix;         /* PREFIX.bam or PREFIX.sam, and PREFIX.fasta */
    dataset_format format;

    size_t         references;
    distribution   ref_length;
    double         covered;    /* fraction of references receiving any reads */
    distribution   reads_per_ref;
    distribution   unmapped;   /* reads aligning nowhere, over the run */

    sim_model      model;

    /* Every value generated derives from this. */
    size_t         seed;
} dataset_config;

/* Writes both files. Returns 0, or -1 with a description in error. */
int dataset_write(const dataset_config *cfg, char *error, size_t error_len);
