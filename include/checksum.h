/* checksum.h -- the digest SAM defines for a reference sequence.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

/* Hex digits in an MD5 digest, which is how an @SQ M5 field spells one. */
#define CHECKSUM_LEN 32

/* Writes the digest of a sequence into hex, which must hold CHECKSUM_LEN + 1 bytes, and
 * terminates it. False only where the digest could not be started.
 *
 * SAM takes the digest over the sequence uppercased and stripped of whitespace. A record as
 * the FASTA reader delivers it holds no whitespace, so case is the only difference, and this
 * is not in general the digest of the bytes as the file holds them. The sequence is
 * uppercased a piece at a time where it lies, a reference running to hundreds of megabytes. */
bool checksum_sequence(const char *seq, size_t len, char *hex);
