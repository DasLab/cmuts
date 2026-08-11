/* shape.c -- the shapes a field's values can have.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "shape.h"

shape_extents shape_per_base(size_t len, size_t cap)
{
    (void)cap;

    return (shape_extents){ .rank = 1, .dim = { len } };
}

shape_extents shape_per_length(size_t len, size_t cap)
{
    (void)len;

    return (shape_extents){ .rank = 1, .dim = { SHAPE_LENGTH_BINS(cap) } };
}

shape_extents shape_none(size_t len, size_t cap)
{
    (void)len;
    (void)cap;

    return (shape_extents){ .rank = 0 };
}

size_t shape_values(shape_fn shape, size_t len, size_t cap)
{
    shape_extents extents = shape(len, cap);
    int           rank    = shape_rank(extents);
    size_t        values  = 1;

    for (int i = 0; i < rank; i++) {
        values *= extents.dim[i];
    }

    return values;
}
