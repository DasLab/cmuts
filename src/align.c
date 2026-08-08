/* align.c -- walking a CIGAR against the reference it was aligned to.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "align.h"

#include "nuc.h"

/* Bits bam_cigar_type sets for an operation that walks along each sequence. */
#define CIGAR_CONSUMES_QUERY     1
#define CIGAR_CONSUMES_REFERENCE 2

static bool is_aligned(uint32_t op)
{
    return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
}

static aln_kind kind_of(uint32_t op)
{
    switch (op) {
        case BAM_CINS:       return ALN_INSERTION;
        case BAM_CDEL:       return ALN_DELETION;
        case BAM_CREF_SKIP:  return ALN_SKIP;
        case BAM_CSOFT_CLIP: return ALN_SOFT_CLIP;
        case BAM_CHARD_CLIP: return ALN_HARD_CLIP;
        default:             return ALN_PAD;
    }
}

/* A position past the end of the reference has nothing to be compared with,
 * which is a thing not known rather than a disagreement. */
static aln_kind compare(const aln_walk *walk, hts_pos_t reference, int32_t query)
{
    nuc theirs;
    nuc ours;

    if (reference < 0 || (size_t)reference >= walk->ref->len)
        return ALN_AMBIGUOUS;

    theirs = nuc_from_char(walk->ref->seq[reference]);
    ours   = nuc_from_read(walk->read->seq, query);

    if (!nuc_is_base(theirs) || !nuc_is_base(ours))
        return ALN_AMBIGUOUS;

    return theirs == ours ? ALN_MATCH : ALN_MISMATCH;
}

/* How far the same answer holds, which is where this run ends and the next
 * begins. */
static uint32_t agreement(const aln_walk *walk, aln_kind kind, uint32_t remaining)
{
    uint32_t len = 1;

    while (len < remaining &&
           compare(walk, walk->reference + (hts_pos_t)len,
                         walk->query + (int32_t)len) == kind)
        len++;

    return len;
}

static void step(aln_walk *walk, uint32_t op, uint32_t len)
{
    int consumes = bam_cigar_type(op);

    walk->consumed += len;

    if (consumes & CIGAR_CONSUMES_QUERY)
        walk->query += (int32_t)len;
    if (consumes & CIGAR_CONSUMES_REFERENCE)
        walk->reference += (hts_pos_t)len;
}

void aln_open(aln_walk *walk, const cm_bam_record *read, const cm_fasta_record *ref)
{
    walk->read      = read;
    walk->ref       = ref;
    walk->op        = 0;
    walk->consumed  = 0;
    walk->reference = read->pos;
    walk->query     = 0;
}

bool aln_next(aln_walk *walk, aln_run *run)
{
    while (walk->op < walk->read->n_cigar) {
        uint32_t cigar     = walk->read->cigar[walk->op];
        uint32_t op        = bam_cigar_op(cigar);
        uint32_t remaining = bam_cigar_oplen(cigar) - walk->consumed;

        if (remaining == 0) {
            walk->op++;
            walk->consumed = 0;
            continue;
        }

        run->reference = walk->reference;
        run->query     = walk->query;

        if (is_aligned(op)) {
            run->kind = compare(walk, walk->reference, walk->query);
            run->len  = agreement(walk, run->kind, remaining);
        } else {
            run->kind = kind_of(op);
            run->len  = remaining;
        }

        step(walk, op, run->len);
        return true;
    }

    return false;
}

/* ------------------------------------------------------------------------ */
/* The path as a band                                                        */
/* ------------------------------------------------------------------------ */

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

    /* A read carrying no CIGAR places nothing, whatever length it stores, and
     * saying so is what keeps a caller from reading centers never written. */
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
