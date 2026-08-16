/* refrow.h -- writing one reference's row from its accumulator.
 *
 * The writer holds the file and the layout; a run holds an accumulator. Two of the fields
 * written are not accumulated: the reactivity and its error are derived here from the mutations
 * and the span, along with the settings governing that.
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

/* Writes every field of one reference's row. Values past len are left at the fill value. */
int refrow_write(refrow *r, int32_t tid, size_t len, const accum *acc,
                 const pairs *pr);
