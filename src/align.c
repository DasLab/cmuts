/* align.c -- where a CIGAR puts each base of a read.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "align.h"

/* Bits bam_cigar_type sets for an operation that walks along each sequence. */
#define CIGAR_CONSUMES_QUERY     1
#define CIGAR_CONSUMES_REFERENCE 2

static bool is_clip(uint32_t op)
{
    return op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP;
}

/* Only a soft clip moves an offset. A hard clipped base is absent from the
 * stored sequence altogether, so the bases either side of it are already
 * neighbours and nothing has to be stepped over. */
static int32_t soft_clipped(uint32_t cigar)
{
    return bam_cigar_op(cigar) == BAM_CSOFT_CLIP
         ? (int32_t)bam_cigar_oplen(cigar)
         : 0;
}

aln_span aln_placed_span(const cm_bam_record *read)
{
    aln_span span  = { .begin = 0, .end = read->l_qseq };
    uint32_t first = 0;
    uint32_t last  = read->n_cigar;

    /* A read carrying no CIGAR places nothing, whatever length it stores.
     * Returning an empty span keeps a caller from reading centers that were
     * never written. */
    if (read->n_cigar == 0)
        return (aln_span){ .begin = 0, .end = 0 };

    while (first < last && is_clip(bam_cigar_op(read->cigar[first])))
        span.begin += soft_clipped(read->cigar[first++]);

    while (last > first && is_clip(bam_cigar_op(read->cigar[last - 1])))
        span.end -= soft_clipped(read->cigar[--last]);

    return span;
}

aln_span aln_places(const cm_bam_record *read, aln_place *places)
{
    aln_span  span      = aln_placed_span(read);
    hts_pos_t reference = read->pos;
    size_t    placed    = 0;

    places[0].first = reference;
    places[0].last  = reference;

    for (uint32_t i = 0; i < read->n_cigar; i++) {
        uint32_t op       = bam_cigar_op(read->cigar[i]);
        uint32_t len      = bam_cigar_oplen(read->cigar[i]);
        int      consumes = bam_cigar_type(op);

        if (is_clip(op))
            continue;

        if (consumes & CIGAR_CONSUMES_QUERY) {
            for (uint32_t j = 0; j < len; j++) {
                if (consumes & CIGAR_CONSUMES_REFERENCE)
                    reference++;

                placed++;
                places[placed].first = reference;
                places[placed].last  = reference;
            }
        } else if (consumes & CIGAR_CONSUMES_REFERENCE) {
            /* Reference passed over with no base to pair it against. The base
             * before it keeps where it was placed and gains where the skip
             * leaves off, so that the stretch it accounts for is the whole of
             * what the path crosses without the read moving. Several skips in a
             * row extend the same one, there being one place per base and not
             * one per operation. */
            reference          += (hts_pos_t)len;
            places[placed].last = reference;
        }
    }

    return span;
}
