/* codec.c -- shuffle and deflate, as HDF5 would have applied them.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "codec.h"

#include <zlib.h>

size_t codec_bound(size_t n_bytes)
{
    return compressBound((uLong)n_bytes);
}

size_t codec_scratch(size_t n_bytes)
{
    return n_bytes;
}

/* Every first byte, then every second, and so on. Runs of values that differ
 * only in their low bytes become runs of equal bytes, which is what makes the
 * deflate that follows worth its cost. */
static void shuffle(unsigned char *dst, const unsigned char *src,
                    size_t n_values, size_t elem_bytes)
{
    for (size_t byte = 0; byte < elem_bytes; byte++) {
        const unsigned char *from = src + byte;

        for (size_t i = 0; i < n_values; i++) {
            *dst++ = *from;
            from  += elem_bytes;
        }
    }
}

size_t codec_encode(unsigned char *dst, unsigned char *scratch,
                    const void *values, size_t n_bytes, size_t elem_bytes)
{
    uLongf written = (uLongf)codec_bound(n_bytes);

    shuffle(scratch, values, n_bytes / elem_bytes, elem_bytes);

    if (compress2(dst, &written, scratch, (uLong)n_bytes,
                  CODEC_DEFLATE_LEVEL) != Z_OK) {
        return 0;
    }

    return written;
}
