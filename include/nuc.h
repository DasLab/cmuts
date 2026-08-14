/* nuc.h -- one nucleotide, decoded from any of the ways it is stored.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stdint.h>

#include <htslib/sam.h>

/* NUC_N is first so that a table indexed by a byte reads as unknown wherever it was not filled
 * in: every IUPAC ambiguity code, every stray character, and the codes BAM uses for a base
 * that could be one of several. */
typedef enum {
    NUC_N,
    NUC_A,
    NUC_C,
    NUC_G,
    NUC_T,
    NUC_COUNT,
} nuc;

/* Named bases, which is how wide a per-base breakdown is. NUC_N is not one of them. */
#define NUC_BASES (NUC_COUNT - 1)

extern const nuc NUC_FROM_CHAR[256];
extern const nuc NUC_FROM_READ[16];

/* Decodes a base of a reference. Case is not significant: lower case in a FASTA marks a
 * repeat, and U reads as T. */
static inline nuc nuc_from_char(char base)
{
    return NUC_FROM_CHAR[(unsigned char)base];
}

/* Decodes a base of a read, by offset into the packed sequence of a cm_bam_record. */
static inline nuc nuc_from_read(const uint8_t *seq, int32_t offset)
{
    return NUC_FROM_READ[bam_seqi(seq, offset)];
}

static inline bool nuc_is_base(nuc base)
{
    return base != NUC_N;
}

/* Gives the index of a named base among NUC_BASES, for addressing a per-base array. */
static inline int nuc_index(nuc base)
{
    return (int)base - (int)NUC_A;
}
