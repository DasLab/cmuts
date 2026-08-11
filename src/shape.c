/* shape.c -- the width and the rank of a row of each kind.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "shape.h"

size_t shape_extent(shape_kind kind, size_t len, size_t cap)
{
    switch (kind) {
        case SHAPE_PER_BASE:   return len;
        case SHAPE_PER_LENGTH: return SHAPE_LENGTH_BINS(cap);
        case SHAPE_SCALAR:     return 1;
        case SHAPE_N_KINDS:    break;
    }

    return 0;
}

int shape_rank(shape_kind kind)
{
    switch (kind) {
        case SHAPE_PER_BASE:
        case SHAPE_PER_LENGTH: return 2;
        case SHAPE_SCALAR:     return 1;
        case SHAPE_N_KINDS:    break;
    }

    return 0;
}
