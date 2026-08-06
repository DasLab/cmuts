/* tally.c -- one read's contribution to a reference.
 *
 * Coverage is real; what counts as a modification is still a placeholder,
 * standing in until that is settled.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "tally.h"

/* Bits bam_cigar_type sets for an operation that walks along each sequence. */
#define CIGAR_CONSUMES_QUERY     1
#define CIGAR_CONSUMES_REFERENCE 2

/* A read and the reference advance at different rates -- an insertion moves
 * along one and not the other -- so neither position follows from the other. */
typedef struct {
    hts_pos_t reference;
    int32_t   query;
} cursor;

typedef struct {
    const cm_bam_record *read;
    const tally_config  *config;
    size_t               limit;   /* bases in the reference */
    accum               *target;
} context;

/* What a base is worth: the chance it was read correctly. A record storing no
 * qualities weighs every base fully, there being nothing to say otherwise. */
static double weight_of(const context *ctx, int32_t query)
{
    return ctx->read->qual
         ? phred_correct(&ctx->config->quality, ctx->read->qual[query])
         : 1.0;
}

/* Coverage counts confidence rather than bases, so a position read poorly is
 * covered less than one read cleanly, and a rate taken against it is divided by
 * the evidence that was really there. */
static void add_coverage(const context *ctx, cursor at, uint32_t len)
{
    double *coverage = accum_data(ctx->target, ACCUM_COVERAGE);

    for (uint32_t i = 0; i < len; i++) {
        hts_pos_t pos = at.reference + (hts_pos_t)i;

        if (pos >= 0 && (size_t)pos < ctx->limit)
            coverage[pos] += weight_of(ctx, at.query + (int32_t)i);
    }
}

static void add_mutations(const context *ctx, hts_pos_t from, uint32_t len)
{
    double *mutations = accum_data(ctx->target, ACCUM_MUTATIONS);

    for (uint32_t i = 0; i < len; i++) {
        hts_pos_t pos = from + (hts_pos_t)i;

        if (pos >= 0 && (size_t)pos < ctx->limit)
            mutations[pos] += 1.0;
    }
}

static void apply_operation(const context *ctx, cursor at, uint32_t op, uint32_t len)
{
    switch (op) {
        case BAM_CMATCH:
        case BAM_CEQUAL:
            add_coverage(ctx, at, len);
            break;
        case BAM_CDIFF:
            add_coverage(ctx, at, len);
            add_mutations(ctx, at.reference, len);
            break;
        case BAM_CDEL:
            add_mutations(ctx, at.reference, len);
            break;
        case BAM_CINS:
            add_mutations(ctx, at.reference, 1);
            break;
        default:
            break;
    }
}

static cursor advance(cursor at, uint32_t op, uint32_t len)
{
    int consumes = bam_cigar_type(op);

    if (consumes & CIGAR_CONSUMES_QUERY)
        at.query += (int32_t)len;
    if (consumes & CIGAR_CONSUMES_REFERENCE)
        at.reference += (hts_pos_t)len;

    return at;
}

void tally_config_build(tally_config *config)
{
    phred_build(&config->quality);
}

void tally(const cm_bam_record *read, const cm_fasta_record *ref,
           const tally_config *config, accum *target)
{
    context ctx = { .read = read, .config = config, .limit = ref->len, .target = target };
    cursor  at  = { .reference = read->pos, .query = 0 };

    for (uint32_t i = 0; i < read->n_cigar; i++) {
        uint32_t op  = bam_cigar_op(read->cigar[i]);
        uint32_t len = bam_cigar_oplen(read->cigar[i]);

        apply_operation(&ctx, at, op, len);
        at = advance(at, op, len);
    }

    *accum_data(target, ACCUM_READS) += 1.0;
}
