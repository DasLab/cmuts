/* metadata.h -- what the output carries besides the per-reference rows.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include "bamstream.h"
#include "h5writer.h"

/* Labels every row with the name of the reference it belongs to, taken from
 * the BAM header rather than from the references actually seen, since a row
 * exists whether or not any read reached it. */
int metadata_write_names(h5writer *out, const cm_bam_stream *bam);

/* Totals belonging to the run as a whole. Unmapped reads are the only such
 * quantity: they align to no reference and so have no row to be counted in.
 * Everything else follows from the per-reference datasets. */
int metadata_write_run(h5writer *out, size_t reads_unmapped);
