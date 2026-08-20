/* chunktally.h -- when an output chunk has received every reference it will.
 *
 * A chunk is finished once the loader has moved past it and every reference it opened
 * there has been written. Either side may be the one that finishes a chunk, so both
 * collect what they finish through chunktally_take_settled.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>
#include <stdint.h>

typedef struct chunktally chunktally;

chunktally *chunktally_create(int32_t n_refs, size_t rows_per_chunk);
void        chunktally_destroy(chunktally *t);

/* Records that the loader opened a reference. References arrive in ascending order,
 * which seals the chunks behind them. */
void chunktally_expect(chunktally *t, int32_t tid);

/* Records that a reference was written. */
void chunktally_wrote(chunktally *t, int32_t tid);

/* Records that nothing more will be opened. */
void chunktally_no_more(chunktally *t);

/* Gives the next finished chunk, or -1 while none is. Each chunk is handed back exactly
 * once. */
int64_t chunktally_take_settled(chunktally *t);
