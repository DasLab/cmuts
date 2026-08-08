/* tally.c -- one read's contribution to a reference.
 *
 * A read whose CIGAR holds an indel is marginalized over the alignments a band
 * around that CIGAR admits, and contributes what each reference position is
 * worth under the posterior rather than under the one path the aligner chose to
 * report. A read without an indel is taken as written: with no gap to move, the
 * posterior sits within rounding of that path at a few hundred times the cost
 * of reading it.
 *
 * Both routes count the same three things. Coverage is the confidence in the
 * bases actually read at a position, so that a rate taken against it is divided
 * by the evidence that was really there; the span is the chance the read
 * reached the position at all, whether it read a base there or passed over it,
 * and is the denominator a deletion belongs over, having no base to be believed.
 * A modification is counted once per event and never once per base: one adduct
 * stops one reverse transcriptase once, whatever length of reference it then
 * skipped, and a disagreement is worth the chance the template really differed
 * rather than the read of it having been wrong. Each event is then worth what
 * its kind is worth, a substitution, a deletion and an insertion not speaking
 * alike of a modification.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "tally.h"

#include <stdlib.h>

#include "align.h"

typedef struct {
    const cm_bam_record   *read;
    const cm_fasta_record *ref;
    const tally_config    *config;
    accum                 *target;
} context;

struct tally_scratch {
    phmm_scratch *phmm;
};

static void add_at(const context *ctx, accum_field_id field, hts_pos_t pos,
                   double value)
{
    double *values = accum_data(ctx->target, field);

    if (pos >= 0 && (size_t)pos < ctx->ref->len)
        values[pos] += value;
}

/* ------------------------------------------------------------------------ */
/* The alignment as written                                                  */
/* ------------------------------------------------------------------------ */

/* A record storing no qualities weighs every base fully, there being nothing to
 * say otherwise. */
static double error_at(const context *ctx, int32_t query)
{
    return ctx->read->qual
         ? phred_error(&ctx->config->quality, ctx->read->qual[query])
         : 0.0;
}

/* Positions the read was set against, whether or not the two could be compared.
 *
 * Coverage counts confidence rather than bases, so a position read poorly is
 * covered less than one read cleanly, and a rate taken against it is divided by
 * the evidence that was really there. The span counts the position whole, being
 * what the read reached rather than what it managed to see there.
 *
 * A comparison that could be made also carries the chance that what was read
 * differed from the reference for a reason other than a misread: almost nothing
 * for a clean agreement, almost one for a clean disagreement. A base neither
 * side has named was still read, and still says nothing either way. */
static void add_aligned_run(const context *ctx, const aln_run *run, aln_kind kind)
{
    const phmm *model = &ctx->config->model;

    for (uint32_t i = 0; i < run->len; i++) {
        hts_pos_t pos   = run->reference + (hts_pos_t)i;
        double    error = error_at(ctx, run->query + (int32_t)i);

        add_at(ctx, ACCUM_COVERAGE, pos, 1.0 - error);
        add_at(ctx, ACCUM_SPANNED, pos, 1.0);

        if (kind != ALN_AMBIGUOUS)
            add_at(ctx, ACCUM_MUTATIONS, pos,
                   phmm_weigh(model, PHMM_SUBSTITUTION,
                              phmm_modification(model, kind == ALN_MATCH,
                                                error)));
    }
}

/* The read reached every position the deletion passed over and read a base at
 * none of them, so it spans them all and covers none. The single event that
 * skipped them is counted where the run ends: reverse transcription reads the
 * template from its 3' end, so the last base passed over is the first the
 * enzyme met. The marginal counts it in the same place. */
static void add_deletion_run(const context *ctx, const aln_run *run)
{
    for (uint32_t i = 0; i < run->len; i++)
        add_at(ctx, ACCUM_SPANNED, run->reference + (hts_pos_t)i, 1.0);

    add_at(ctx, ACCUM_MUTATIONS,
           run->reference + (hts_pos_t)run->len - 1,
           phmm_weigh(&ctx->config->model, PHMM_DELETION, 1.0));
}

/* An insertion sits between two reference positions rather than on one, so it
 * is counted once at the position it precedes however many bases it carries,
 * and reaches neither. */
static void add_insertion_run(const context *ctx, const aln_run *run)
{
    add_at(ctx, ACCUM_MUTATIONS, run->reference,
           phmm_weigh(&ctx->config->model, PHMM_INSERTION, 1.0));
}

/* A clip, a skip and a pad reach no reference position and are counted nowhere;
 * they are reported so that the walk stays in step, not so that it accumulates.
 */
static void apply_run(const context *ctx, const aln_run *run)
{
    switch (run->kind) {
        case ALN_MATCH:
        case ALN_MISMATCH:
        case ALN_AMBIGUOUS: add_aligned_run(ctx, run, run->kind); break;
        case ALN_DELETION:  add_deletion_run(ctx, run);           break;
        case ALN_INSERTION: add_insertion_run(ctx, run);          break;
        default:                                                  break;
    }
}

static void walk_alignment(const context *ctx)
{
    aln_walk walk;
    aln_run  run;

    aln_open(&walk, ctx->read, ctx->ref);

    while (aln_next(&walk, &run))
        apply_run(ctx, &run);
}

/* ------------------------------------------------------------------------ */
/* The alignment marginalized                                                */
/* ------------------------------------------------------------------------ */

static bool has_indel(const cm_bam_record *read)
{
    for (uint32_t i = 0; i < read->n_cigar; i++) {
        uint32_t op = bam_cigar_op(read->cigar[i]);

        if (op == BAM_CINS || op == BAM_CDEL)
            return true;
    }

    return false;
}

static void add_window(const context *ctx, const phmm_window *window)
{
    for (size_t i = 0; i < window->len; i++) {
        hts_pos_t pos = window->origin + (hts_pos_t)i;

        add_at(ctx, ACCUM_COVERAGE, pos, window->coverage[i]);
        add_at(ctx, ACCUM_SPANNED, pos, window->spanned[i]);
        add_at(ctx, ACCUM_MUTATIONS, pos, window->mutations[i]);
    }
}

/* Only a read the aligner had a choice about is worth marginalizing, and only
 * an indel gives it one. */
static bool marginalize(const context *ctx, tally_scratch *scratch)
{
    phmm_window window;

    if (!scratch || !has_indel(ctx->read))
        return false;

    if (!phmm_run(&ctx->config->model, &ctx->config->quality,
                  ctx->read, ctx->ref, scratch->phmm, &window))
        return false;

    add_window(ctx, &window);

    return true;
}

/* ------------------------------------------------------------------------ */
/* Setup                                                                     */
/* ------------------------------------------------------------------------ */

void tally_config_build(tally_config *config)
{
    phmm_params  params  = phmm_defaults();
    phmm_weights weights = phmm_default_weights();

    phred_build(&config->quality);
    phmm_build(&config->model, &params, &weights);
}

tally_scratch *tally_scratch_create(void)
{
    tally_scratch *scratch = calloc(1, sizeof *scratch);

    if (!scratch)
        return NULL;

    scratch->phmm = phmm_scratch_create();

    if (!scratch->phmm) {
        free(scratch);
        return NULL;
    }

    return scratch;
}

void tally_scratch_destroy(tally_scratch *scratch)
{
    if (!scratch)
        return;

    phmm_scratch_destroy(scratch->phmm);
    free(scratch);
}

void tally(const cm_bam_record *read, const cm_fasta_record *ref,
           const tally_config *config, tally_scratch *scratch, accum *target)
{
    context ctx = {
        .read   = read,
        .ref    = ref,
        .config = config,
        .target = target,
    };

    if (!marginalize(&ctx, scratch))
        walk_alignment(&ctx);

    *accum_data(target, ACCUM_READS) += 1.0;
}
