/* shape.h -- the extents one reference's values occupy.
 *
 * A field is described by a function from the reference it belongs to, and the longest in
 * the run, to the extents of the values it holds there. The rank is how many extents there
 * are, so it is answered by the same call rather than by a second description that could
 * disagree with the first. A field of one value per reference has no extents at all: the
 * empty vector, whose product is the single value it holds.
 *
 * The accumulator and the output both describe their fields this way, and must, since a
 * value accumulated per read length has to reach a dataset column of the same meaning. What
 * the output adds is the reference dimension, which an accumulator holding one reference
 * has no use for.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

/* Bins a read-length histogram covers, given the longest reference in the run: one for
 * every length from 1 to twice it, so that bin i holds the reads of length i + 1. The
 * range reaches past a reference because a read carrying insertions or soft-clipped ends
 * is longer than the one it aligns to. It starts at 1 because a read storing no sequence
 * is refused before it is counted, leaving no read of length zero to hold.
 *
 * Every row is this wide, whatever its own reference measures. A read length is not a
 * position in a reference, so a column a short reference has no reads for is a count of
 * zero and not padding; sizing rows individually would leave the same column meaning
 * different things in different rows.
 *
 * A read longer than the range is counted in no bin. How many there were is the reads
 * total less the row's own sum. */
#define SHAPE_LENGTH_BINS(cap) (2 * (cap))

/* The largest number of extents any of the shapes below writes. Raise it alongside a shape
 * that writes more. */
#define SHAPE_RANK_MAX 1

/* The extents one reference's values occupy. A shape fills this and nothing wider, so one
 * reporting a rank it has no room for misdescribes a field rather than overrunning
 * anything; shape_rank is what every reader of dim counts through. */
typedef struct {
    int    rank;
    size_t dim[SHAPE_RANK_MAX];
} shape_extents;

/* Gives the extents a reference of len bases occupies, in a run whose longest is cap.
 *
 * The rank is the same whatever a shape is given; only the extents vary, which is what
 * lets the rank of a dataset be asked for without a run to ask it of. */
typedef shape_extents (*shape_fn)(size_t len, size_t cap);

shape_extents shape_per_base(size_t len, size_t cap);
shape_extents shape_per_length(size_t len, size_t cap);
shape_extents shape_none(size_t len, size_t cap);

/* What a shape is called, for describing a field to something outside the program. A shape
 * added above is named here as well, or it describes itself as unknown. */
const char *shape_name(shape_fn shape);

/* Gives the extents held: what was reported, or all there is room for, whichever is fewer.
 * Inline so that every loop over dim carries the bound with it. */
static inline int shape_rank(shape_extents extents)
{
    return extents.rank < SHAPE_RANK_MAX ? extents.rank : SHAPE_RANK_MAX;
}

/* The values one reference occupies, the product of its extents. */
size_t shape_values(shape_fn shape, size_t len, size_t cap);
