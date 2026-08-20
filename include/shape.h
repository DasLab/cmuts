/* shape.h -- the extents one reference's values occupy.
 *
 * A shape maps a reference's length, and the longest in the run, to the extents the
 * field's values occupy there. The accumulator and the output describe their fields
 * with the same shapes.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

/* The bins a read-length histogram holds: one for each length from 1 to twice the
 * longest reference. */
#define SHAPE_LENGTH_BINS(cap) (2 * (cap))

/* The largest number of extents any of the shapes below writes. Raise it alongside a shape
 * that writes more. */
#define SHAPE_RANK_MAX 2

/* The extents one reference's values occupy. Loops over dim take their bound from
 * shape_rank. */
typedef struct {
    int    rank;
    size_t dim[SHAPE_RANK_MAX];
} shape_extents;

/* Gives the extents a reference of len bases occupies, in a run whose longest is cap.
 * The rank does not depend on the arguments, so a call at any length gives a dataset's
 * rank. */
typedef shape_extents (*shape_fn)(size_t len, size_t cap);

shape_extents shape_per_base(size_t len, size_t cap);
shape_extents shape_per_pair(size_t len, size_t cap);
shape_extents shape_per_length(size_t len, size_t cap);
shape_extents shape_none(size_t len, size_t cap);

/* Gives the name of a shape, for the JSON description of a field, or "unknown" for a
 * shape without one. */
const char *shape_name(shape_fn shape);

/* Gives extent i of a shape as the documentation writes it, "l" standing for the longest
 * reference, or NULL where the extent does not vary with the run. */
const char *shape_symbol(shape_fn shape, int extent);

/* Gives the rank, capped at SHAPE_RANK_MAX. */
static inline int shape_rank(shape_extents extents)
{
    return extents.rank < SHAPE_RANK_MAX ? extents.rank : SHAPE_RANK_MAX;
}

/* Gives the values one reference occupies, which is the product of its extents. */
size_t shape_values(shape_fn shape, size_t len, size_t cap);
