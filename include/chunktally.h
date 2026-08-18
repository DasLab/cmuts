/* chunktally.h -- when an output chunk has received every reference it will.
 *
 * A chunk cannot be written until nothing further can land in it, and how long
 * that takes is not bounded: a reference whose reads are slow to process
 * finishes arbitrarily far behind ones opened after it. So the question is
 * answered by counting rather than by waiting a fixed distance -- a chunk is
 * finished once the loader has moved past it and every reference it opened
 * there has been written.
 *
 * The loader announces; the consumer reports; either may be the one that
 * settles a chunk, so both hand back what they settle through take_settled.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdint.h>

typedef struct chunktally chunktally;

chunktally *chunktally_create(int32_t n_refs, size_t rows_per_chunk);
void        chunktally_destroy(chunktally *t);

/* The loader, as each reference is opened. They arrive in ascending order,
 * which is what seals the chunks behind them. */
void chunktally_expect(chunktally *t, int32_t tid);

/* The consumer, as each reference is written. */
void chunktally_wrote(chunktally *t, int32_t tid);

/* The loader, once the file is spent: nothing more will be opened. */
void chunktally_no_more(chunktally *t);

/* The next chunk that has become finished, or -1 while none has. Each chunk is
 * handed back exactly once. */
int64_t chunktally_take_settled(chunktally *t);
