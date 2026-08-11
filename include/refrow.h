/* refrow.h -- writing one reference's row from its accumulator.
 *
 * The writer holds the file and the layout; a run holds an accumulator. Two of
 * the fields written are not accumulated: the reactivity and its error are
 * derived from the mutations and the span, which is done here, along with the
 * settings that govern it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdint.h>

#include "accum.h"
#include "h5writer.h"
#include "rates.h"

typedef struct refrow refrow;

/* Borrows the writer, which must outlive it. */
refrow *refrow_create(h5writer *out, rate_config rates, size_t ref_cap);
void    refrow_destroy(refrow *r);

/* Writes every field of one reference's row. Values past len are left at the
 * fill value. */
int refrow_write(refrow *r, int32_t tid, size_t len, const accum *acc);
