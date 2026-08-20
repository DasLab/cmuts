/* checksum.c -- MD5 over a reference sequence, as SAM defines it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "checksum.h"

#include <ctype.h>

#include <htslib/hts.h>

/* Bases uppercased at a time: enough to feed the digest in few pieces, small enough to
 * sit on the stack. */
#define CHUNK_BASES 4096

static void hash_uppercase(hts_md5_context *md5, const char *seq, size_t len,
                           char from, char to)
{
    char chunk[CHUNK_BASES];

    for (size_t done = 0; done < len; ) {
        size_t n = len - done < CHUNK_BASES ? len - done : CHUNK_BASES;

        for (size_t i = 0; i < n; i++) {
            char base = (char)toupper((unsigned char)seq[done + i]);

            chunk[i] = base == from ? to : base;
        }

        hts_md5_update(md5, chunk, n);
        done += n;
    }
}

bool checksum_sequence(const char *seq, size_t len, char *hex)
{
    return checksum_sequence_swapped(seq, len, 0, 0, hex);
}

bool checksum_sequence_swapped(const char *seq, size_t len, char from, char to,
                               char *hex)
{
    hts_md5_context *md5 = hts_md5_init();
    unsigned char    digest[CHECKSUM_LEN / 2];

    if (!md5) {
        return false;
    }

    hash_uppercase(md5, seq, len, from, to);
    hts_md5_final(digest, md5);
    hts_md5_hex(hex, digest);
    hts_md5_destroy(md5);

    return true;
}
