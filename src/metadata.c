/* metadata.c -- reference names and run totals.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "metadata.h"

#include <stdlib.h>

#define ATTRIBUTE_READS_UNMAPPED "reads_unmapped"

int metadata_write_names(h5writer *out, const cm_bam_reader *bam)
{
    int32_t      n     = cm_bam_nref(bam);
    const char **names = calloc((size_t)n, sizeof *names);
    int          status;

    if (!names)
        return -1;

    /* Borrowed from the header, which outlives the write. */
    for (int32_t tid = 0; tid < n; tid++)
        names[tid] = cm_bam_refname(bam, tid);

    status = h5writer_names(out, names, n);
    free(names);

    return status;
}

int metadata_write_run(h5writer *out, size_t reads_unmapped)
{
    return h5writer_count(out, ATTRIBUTE_READS_UNMAPPED, reads_unmapped);
}
