/* tally.c -- one read's contribution to a reference.
 *
 * Every read is marginalized over the alignments a band around its CIGAR
 * admits, and contributes what each reference position is worth under the
 * posterior rather than under the one path the aligner chose to report. A band
 * of nothing is that path and nothing else, so taking an alignment as written
 * is the narrowest setting rather than a second way of counting.
 *
 * Three things are counted. Coverage is the confidence in the bases actually
 * read at a position, so that a rate taken against it is divided by the
 * evidence that was really there; the span is the chance the read reached the
 * position at all, whether it read a base there or passed over it, and is the
 * denominator a deletion belongs over, having no base to be believed. A
 * modification is counted once per event and never once per base: one adduct
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

typedef struct {
    const cm_bam_record   *read;
    const cm_fasta_record *ref;
    const tally_tables    *tables;
    accum                 *target;
} context;

struct tally_scratch {
    phmm_scratch *phmm;
    int          *half;   /* the band handed to the marginal, row by row */
    size_t        rows;   /* rows it is sized for */
};

static void add_at(const context *ctx, accum_field_id field, hts_pos_t pos,
                   double value)
{
    double *values = accum_data(ctx->target, field);

    if (pos >= 0 && (size_t)pos < ctx->ref->len)
        values[pos] += value;
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

/* One band, the same width at every row, which is the shape the marginal was
 * built around and the only one asked for yet. It is grown to the longest read
 * seen and filled once, the width being settled before any read arrives.
 *
 * Returns NULL where it cannot be grown, which is the run over: there is no
 * other way a read is counted. */
static const int *uniform_band(tally_scratch *scratch, const cm_bam_record *read,
                               int band)
{
    size_t rows = (size_t)read->l_qseq + 1;
    int   *half;

    if (rows <= scratch->rows)
        return scratch->half;

    half = realloc(scratch->half, rows * sizeof *half);

    if (!half)
        return NULL;

    for (size_t i = 0; i < rows; i++)
        half[i] = band;

    scratch->half = half;
    scratch->rows = rows;

    return half;
}

static phmm_status marginalize(const context *ctx, tally_scratch *scratch)
{
    phmm_window window;
    phmm_status status;
    const int  *half = uniform_band(scratch, ctx->read, ctx->tables->band);

    if (!half)
        return PHMM_NO_MEMORY;

    status = phmm_run(&ctx->tables->model, &ctx->tables->quality,
                      ctx->read, ctx->ref, half, scratch->phmm, &window);

    if (status == PHMM_OK)
        add_window(ctx, &window);

    return status;
}

/* ------------------------------------------------------------------------ */
/* Setup                                                                     */
/* ------------------------------------------------------------------------ */

tally_config tally_defaults(void)
{
    return (tally_config){ .band = PHMM_DEFAULT_BAND };
}

void tally_tables_build(tally_tables *tables, const tally_config *config)
{
    phmm_params  params  = phmm_defaults();
    phmm_weights weights = phmm_default_weights();

    phred_build(&tables->quality);
    phmm_build(&tables->model, &params, &weights);
    tables->band = config->band;
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
    free(scratch->half);
    free(scratch);
}

phmm_status tally(const cm_bam_record *read, const cm_fasta_record *ref,
                  const tally_tables *tables, tally_scratch *scratch,
                  accum *target)
{
    context ctx = {
        .read   = read,
        .ref    = ref,
        .tables = tables,
        .target = target,
    };
    phmm_status status = marginalize(&ctx, scratch);

    if (status != PHMM_OK)
        return status;

    *accum_data(target, ACCUM_READS) += 1.0;

    return PHMM_OK;
}
