/* codec.h -- the filter chain a chunk carries, applied by hand.
 *
 * A chunk handed to HDF5 whole arrives already filtered, so the shuffle and
 * deflate the library would otherwise apply are applied here instead, in the
 * order h5layout_creation_plist declares them. That declaration and this file
 * are two halves of one contract: a dataset filtered any other way cannot be
 * written through here. The deflate level lives here and h5layout passes it to
 * H5Pset_deflate, so the two halves cannot drift apart.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stddef.h>

#define CODEC_DEFLATE_LEVEL 3

/* The most codec_encode can write for values of these bytes. */
size_t codec_bound(size_t n_bytes);

/* Scratch of this size lets one thread encode without allocating. */
size_t codec_scratch(size_t n_bytes);

/* Filters values into dst, which holds codec_bound(n_bytes), by way of
 * scratch, which holds codec_scratch(n_bytes). elem_bytes is the width of one
 * value, which the shuffle works in. Returns the bytes written, or zero if the
 * deflate failed. */
size_t codec_encode(unsigned char *dst, unsigned char *scratch,
                    const void *values, size_t n_bytes, size_t elem_bytes);
