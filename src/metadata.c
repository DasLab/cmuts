/* metadata.c -- run totals.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "metadata.h"

/* Named within the group, which supplies the reads the name drops. */
#define DATASET_READS_UNMAPPED "unmapped"

int metadata_write_run(h5writer *out, size_t reads_unmapped)
{
    return h5writer_count(out, DATASET_READS_UNMAPPED, reads_unmapped);
}

int metadata_read_run(h5reader *in, size_t *reads_unmapped)
{
    return h5reader_count(in, DATASET_READS_UNMAPPED, reads_unmapped);
}
