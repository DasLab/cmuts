/* tally.c -- one read's contribution to a reference.
 *
 * Coverage is real. What counts as a modification is still a placeholder: a
 * disagreement of any kind counts once per reference position it covers, and
 * an insertion once at the position it sits before.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "tally.h"

#include "align.h"

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
static void add_coverage(const context *ctx, const aln_run *run)
{
    double *coverage = accum_data(ctx->target, ACCUM_COVERAGE);

    for (uint32_t i = 0; i < run->len; i++) {
        hts_pos_t pos = run->reference + (hts_pos_t)i;

        if (pos >= 0 && (size_t)pos < ctx->limit)
            coverage[pos] += weight_of(ctx, run->query + (int32_t)i);
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

/* An insertion sits between two reference positions rather than on one, so it
 * is recorded once at the position it precedes however many bases it carries.
 * A clip, a skip and a pad cover no reference position and are counted nowhere;
 * they are reported so that the walk stays in step, not so that it accumulates.
 */
static void apply_run(const context *ctx, const aln_run *run)
{
    switch (run->kind) {
        case ALN_MATCH:
        case ALN_AMBIGUOUS:
            add_coverage(ctx, run);
            break;
        case ALN_MISMATCH:
            add_coverage(ctx, run);
            add_mutations(ctx, run->reference, run->len);
            break;
        case ALN_DELETION:
            add_mutations(ctx, run->reference, run->len);
            break;
        case ALN_INSERTION:
            add_mutations(ctx, run->reference, 1);
            break;
        default:
            break;
    }
}

void tally_config_build(tally_config *config)
{
    phred_build(&config->quality);
}

void tally(const cm_bam_record *read, const cm_fasta_record *ref,
           const tally_config *config, accum *target)
{
    context  ctx = { .read = read, .config = config, .limit = ref->len, .target = target };
    aln_walk walk;
    aln_run  run;

    aln_open(&walk, read, ref);

    while (aln_next(&walk, &run))
        apply_run(&ctx, &run);

    *accum_data(target, ACCUM_READS) += 1.0;
}
