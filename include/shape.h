/* shape.h -- how many values a per-reference field occupies.
 *
 * The accumulator and the output both describe their fields by the kind of row each has,
 * and the two must agree on how wide a row of a given kind is: a value accumulated per
 * read length has to land in a dataset column of the same meaning. The kinds and their
 * widths are stated here once, and neither module states them again.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

typedef enum {
    SHAPE_PER_BASE,    /* one value per reference base */
    SHAPE_PER_LENGTH,  /* one value per read length, plus an overflow bin */
    SHAPE_SCALAR,      /* one value per reference */
    SHAPE_N_KINDS,
} shape_kind;

/* Bins a read-length histogram covers, given the longest reference in the run: one for
 * every length from 0 to twice it. The range reaches past a reference because a read
 * carrying insertions or soft-clipped ends is longer than the one it aligns to.
 *
 * Every row is this wide, whatever its own reference measures. A read length is not a
 * position in a reference, so a column a short reference has no reads for is a count of
 * zero and not padding; sizing rows individually would leave the same column meaning
 * different things in different rows.
 *
 * A read longer than the range is counted in no bin. How many there were is the reads
 * total less the row's own sum. */
#define SHAPE_LENGTH_BINS(cap) (2 * (cap) + 1)

/* Values a field of this kind occupies for a reference of len bases, in a run whose
 * longest is cap. */
size_t shape_extent(shape_kind kind, size_t len, size_t cap);

/* Dimensions a dataset of this kind has: the reference it is indexed by, and whatever the
 * kind adds to that. */
int shape_rank(shape_kind kind);

/* The largest rank any kind has, which is how long an array describing a dataset of any of
 * them must be. Raise it alongside a kind of a rank not yet seen. */
#define SHAPE_RANK_MAX 2
