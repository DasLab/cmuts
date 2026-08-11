/* metadata.h -- what the output carries besides the per-reference rows.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#include "h5reader.h"
#include "h5writer.h"

/* Totals belonging to the run as a whole. Unmapped reads are the only such
 * quantity: they align to no reference and so have no row to be counted in.
 * Everything else follows from the per-reference datasets. */
int metadata_write_run(h5writer *out, size_t reads_unmapped);

/* The same totals, read back from a file already written. */
int metadata_read_run(h5reader *in, size_t *reads_unmapped);
