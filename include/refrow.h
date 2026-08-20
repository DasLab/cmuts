/* refrow.h -- writing one reference's row from its accumulator.
 *
 * The writer holds the file and the layout; a run holds an accumulator. The fields not
 * accumulated directly are derived here: the reactivity and its error from the mutations
 * and the span, and the pairwise statistics from the pair sums.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdint.h>

#include "accum.h"
#include "h5writer.h"
#include "pairs.h"
#include "rates.h"

typedef struct refrow refrow;

/* Borrows the writer, which must outlive it. wanted names the optional fields the run
 * writes, and sizes the scratch a derived row is computed into. */
refrow *refrow_create(h5writer *out, rate_config rates, size_t ref_cap,
                      const bool *wanted);
void    refrow_destroy(refrow *r);

/* Writes one reference's row: for every field, its values where the reference has any,
 * and the region outside the reference where the field's fill does not already read as
 * outside it. Anything else keeps the fill.
 *
 * acc is NULL for a reference no read arrived on, which has no accumulated values. seq is
 * its sequence, which it has whether or not a read arrived. */
int refrow_write(refrow *r, int32_t tid, size_t len, const char *seq, const accum *acc,
                 const pairs *pr);
