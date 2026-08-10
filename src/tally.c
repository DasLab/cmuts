/* tally.c -- one read's contribution to a reference.
 *
 * Drives the marginal over one read and adds the window it returns to the
 * reference's accumulator, clipping to the reference's own bounds. What the
 * three counted quantities mean is described in phmm.c, which computes them,
 * and in accum.c, which stores them.
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

/* The read's own stored length, which is what the length filters are applied
 * to and so counts inserted and soft-clipped bases. A read longer than the
 * range the bins cover is counted in none of them; the reads total is what
 * says how many those were. */
static void add_length(const context *ctx)
{
    double *bins   = accum_data(ctx->target, ACCUM_LENGTHS);
    size_t  length = (size_t)ctx->read->l_qseq;

    if (length < ACCUM_LENGTH_BINS(ctx->ref->len))
        bins[length] += 1.0;
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

/* A band of uniform width, the only shape used so far. Grown to the
 * longest read seen and filled once, the width being fixed before any read
 * arrives. Returns NULL if it cannot be grown, which ends the run. */
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
    add_length(&ctx);

    return PHMM_OK;
}
