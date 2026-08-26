/* tally.c -- one read's contribution to a reference.
 *
 * Runs the marginal over one read and adds the window it returns to the reference's
 * accumulator, clipping to the reference's bounds. The counted quantities are
 * described in phmm.h.
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
    pairs                 *target_pairs;   /* NULL where no pairs are counted */
} context;

struct tally_scratch {
    phmm_scratch *phmm;
    int          *half;   /* the band passed to the marginal, one per row */
    size_t        rows;   /* rows it is sized for */
};

/* Counts the read in its length bin. Binned by stored length, as the length filters are
 * applied, so inserted and soft-clipped bases count. A read longer than the range the
 * bins cover falls in none of them; the reads total gives how many those were.
 *
 * The bins begin at length 1, since a read storing no sequence has been refused
 * already, so the guard against zero is against that filter changing and not against
 * anything reachable from here. */
static void add_length(const context *ctx)
{
    double *bins   = accum_data(ctx->target, ACCUM_LENGTHS);
    size_t  length = (size_t)ctx->read->l_qseq;

    if (length > 0 && length <= SHAPE_LENGTH_BINS(ctx->target->cap)) {
        bins[length - 1] += 1.0;
    }
}

/* Adds the window's stretch on the reference to the per-base fields, clipping once
 * rather than testing each position. The window holds one read, so each event value is
 * the posterior chance this read carries an event of that kind at the position. */
static void add_window(const context *ctx, const phmm_window *window)
{
    double *coverage   = accum_data(ctx->target, ACCUM_COVERAGE);
    double *mismatches = accum_data(ctx->target, ACCUM_MISMATCHES);
    double *insertions = accum_data(ctx->target, ACCUM_INSERTIONS);
    double *deletions  = accum_data(ctx->target, ACCUM_DELETIONS);
    double *ends       = accum_data(ctx->target, ACCUM_ENDS);
    size_t  begin;
    size_t  end;

    phmm_window_bounds(window, ctx->ref->len, &begin, &end);

    for (size_t i = begin; i < end; i++) {
        size_t pos = (size_t)(window->origin + (hts_pos_t)i);

        coverage[pos]   += window->coverage[i];
        mismatches[pos] += window->mismatches[i];
        insertions[pos] += window->insertions[i];
        deletions[pos]  += window->deletions[i];
        ends[pos]       += window->ends[i];
    }
}

/* Gives a band of uniform width, the only shape used so far. It is grown to the longest
 * read seen and filled once, the width being fixed before any read arrives. Returns
 * NULL if it cannot be grown, which ends the run. */
static const int *uniform_band(tally_scratch *scratch, const cm_bam_record *read,
                               int band)
{
    size_t rows = (size_t)read->l_qseq + 1;
    int   *half;

    if (rows <= scratch->rows) {
        return scratch->half;
    }

    half = realloc(scratch->half, rows * sizeof *half);

    if (!half) {
        return NULL;
    }

    for (size_t i = 0; i < rows; i++) {
        half[i] = band;
    }

    scratch->half = half;
    scratch->rows = rows;

    return half;
}

static phmm_status marginalize(const context *ctx, tally_scratch *scratch)
{
    phmm_window window;
    phmm_status status;
    const int  *half = uniform_band(scratch, ctx->read, ctx->tables->band);

    if (!half) {
        return PHMM_NO_MEMORY;
    }

    status = phmm_run(&ctx->tables->model, &ctx->tables->profile,
                      &ctx->tables->quality, ctx->read, ctx->ref, half,
                      scratch->phmm, &window);

    if (status == PHMM_OK) {
        add_window(ctx, &window);

        if (ctx->target_pairs) {
            pairs_count(ctx->target_pairs, ctx->ref->len, &window);
        }
    }

    return status;
}

/* ------------------------------------------------------------------------ */
/* Setup                                                                     */
/* ------------------------------------------------------------------------ */

tally_config tally_defaults(void)
{
    return (tally_config){
        .band      = PHMM_DEFAULT_BAND,
        .min_phred = 0,
        .params    = phmm_defaults(),
    };
}

/* Fills one profile array with a single rate. */
static void fill_uniform(double *values, size_t cap, double rate)
{
    for (size_t i = 0; i < cap; i++) {
        values[i] = rate;
    }
}

int tally_tables_build(tally_tables *tables, const tally_config *config, size_t cap)
{
    const phmm_params *params = &config->params;
    double            *rates  = malloc(3 * cap * sizeof *rates);

    if (!rates) {
        return -1;
    }

    phred_build(&tables->quality, config->min_phred);
    phmm_build(&tables->model, params);
    tables->band  = config->band;
    tables->rates = rates;

    tables->profile.modification   = rates;
    tables->profile.open_insertion = rates + cap;
    tables->profile.open_deletion  = rates + 2 * cap;

    fill_uniform(rates, cap, params->modification);
    fill_uniform(rates + cap, cap, params->open_insertion);
    fill_uniform(rates + 2 * cap, cap, params->open_deletion);

    return 0;
}

void tally_tables_free(tally_tables *tables)
{
    free(tables->rates);
    tables->rates = NULL;
}

tally_scratch *tally_scratch_create(void)
{
    tally_scratch *scratch = calloc(1, sizeof *scratch);

    if (!scratch) {
        return NULL;
    }

    scratch->phmm = phmm_scratch_create();

    if (!scratch->phmm) {
        free(scratch);
        return NULL;
    }

    return scratch;
}

void tally_scratch_destroy(tally_scratch *scratch)
{
    if (!scratch) {
        return;
    }

    phmm_scratch_destroy(scratch->phmm);
    free(scratch->half);
    free(scratch);
}

phmm_status tally(const cm_bam_record *read, const cm_fasta_record *ref,
                  const tally_tables *tables, tally_scratch *scratch,
                  accum *target, pairs *target_pairs)
{
    context ctx = {
        .read         = read,
        .ref          = ref,
        .tables       = tables,
        .target       = target,
        .target_pairs = target_pairs,
    };
    phmm_status status = marginalize(&ctx, scratch);

    /* A read the model gives no path is counted where a read a filter turned away is
     * counted, and adds no other value. No value has reached the target: the window is added
     * only on PHMM_OK, and the failure is seen before either count below. */
    if (status == PHMM_NO_PATH) {
        *accum_data(target, ACCUM_FILTERED) += 1.0;
        return status;
    }

    if (status != PHMM_OK) {
        return status;
    }

    *accum_data(target, ACCUM_READS) += 1.0;
    add_length(&ctx);

    return PHMM_OK;
}
