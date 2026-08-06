/* process.c -- placeholder read processing.
 *
 * Stands in for the real analysis while the load-sharing machinery is built
 * out. It walks the CIGAR and counts aligned and non-matching reference
 * positions, which is enough to give every reference a non-trivial, exactly
 * reproducible result to check the parallel merge against.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "process.h"

/* Bit set by bam_cigar_type for operations that advance along the reference. */
#define CIGAR_CONSUMES_REFERENCE 2

/* Arithmetic burned per read, standing in for the cost of real processing.
 * Zero for ordinary runs; raise it to put the worker pool under load and
 * measure how it scales. Goes away with the placeholder. */
#define SYNTHETIC_LOAD 0

/* The result is deposited in a local volatile so the compiler cannot discard
 * the loop, and depends only on the seed, so it stays deterministic and
 * thread-independent: the accumulators never see it. */
static void burn(hts_pos_t seed)
{
    volatile double sink;
    double          x = (double)(seed | 1);

    for (unsigned i = 0; i < SYNTHETIC_LOAD; i++)
        x = x * 1.0000001 + 1.0;

    sink = x;
    (void)sink;
}

static void add_span(double *field, hts_pos_t from, uint32_t len, size_t limit)
{
    for (uint32_t i = 0; i < len; i++) {
        hts_pos_t pos = from + (hts_pos_t)i;

        if (pos >= 0 && (size_t)pos < limit)
            field[pos] += 1.0;
    }
}

static void apply_operation(uint32_t op, uint32_t len, hts_pos_t pos,
                            const cm_fasta_record *ref, accum *target)
{
    double *coverage  = accum_data(target, ACCUM_COVERAGE);
    double *mutations = accum_data(target, ACCUM_MUTATIONS);

    switch (op) {
        case BAM_CMATCH:
        case BAM_CEQUAL:
            add_span(coverage, pos, len, ref->len);
            break;
        case BAM_CDIFF:
            add_span(coverage, pos, len, ref->len);
            add_span(mutations, pos, len, ref->len);
            break;
        case BAM_CDEL:
            add_span(mutations, pos, len, ref->len);
            break;
        case BAM_CINS:
            add_span(mutations, pos, 1, ref->len);
            break;
        default:
            break;
    }
}

void process(const cm_bam_record *read, const cm_fasta_record *ref, accum *target)
{
    hts_pos_t pos = read->pos;

    for (uint32_t i = 0; i < read->n_cigar; i++) {
        uint32_t op  = bam_cigar_op(read->cigar[i]);
        uint32_t len = bam_cigar_oplen(read->cigar[i]);

        apply_operation(op, len, pos, ref, target);

        if (bam_cigar_type(op) & CIGAR_CONSUMES_REFERENCE)
            pos += len;
    }

    *accum_data(target, ACCUM_READS) += 1.0;

    burn(read->pos);
}
