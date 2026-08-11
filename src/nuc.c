/* nuc.c -- the tables the two encodings are read through.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "nuc.h"

const nuc NUC_FROM_CHAR[256] = {
    ['A'] = NUC_A, ['a'] = NUC_A,
    ['C'] = NUC_C, ['c'] = NUC_C,
    ['G'] = NUC_G, ['g'] = NUC_G,
    ['T'] = NUC_T, ['t'] = NUC_T,
    ['U'] = NUC_T, ['u'] = NUC_T,
};

/* BAM gives each base four bits, one per possibility, so only the four powers of two
 * name a single base. */
const nuc NUC_FROM_READ[16] = {
    [1] = NUC_A,
    [2] = NUC_C,
    [4] = NUC_G,
    [8] = NUC_T,
};

char nuc_to_char(nuc base)
{
    static const char LETTER[NUC_COUNT] = {
        [NUC_N] = 'N',
        [NUC_A] = 'A',
        [NUC_C] = 'C',
        [NUC_G] = 'G',
        [NUC_T] = 'T',
    };

    return base < NUC_COUNT ? LETTER[base] : 'N';
}
