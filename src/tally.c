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

#include "filter.h"

/* A value of ends equal at every base, which leaves every placement of a read's 5'-most
 * paired base within its band equally likely. */
#define UNIFORM_ENDS 1.0

typedef struct {
    const cm_bam_record   *read;
    const cm_fasta_record *ref;
    const phmm_rates      *rates;
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
 * A supplementary alignment stores the piece it places and not the read, so binning one
 * would record a length no molecule had.
 *
 * The bins begin at length 1, since a read storing no sequence has been refused
 * already, so the guard against zero is against that filter changing and not against
 * anything reachable from here. */
static void add_length(const context *ctx)
{
    double *bins;
    size_t  length;

    if (filter_is_supplementary(ctx->read)) {
        return;
    }

    bins   = accum_data(ctx->target, ACCUM_LENGTHS);
    length = (size_t)ctx->read->l_qseq;

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
        size_t               pos = (size_t)(window->origin + (hts_pos_t)i);
        const phmm_position *at  = &window->at[i];

        coverage[pos]   += at->coverage;
        mismatches[pos] += at->mismatches;
        insertions[pos] += at->insertions;
        deletions[pos]  += at->deletions;
        ends[pos]       += at->ends;
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

    status = phmm_run(ctx->rates, &ctx->tables->quality, ctx->read, ctx->ref, half,
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
        .uniform   = phmm_uniform_defaults(),
    };
}

/* Fills one array of rates with a single rate. */
static void fill_uniform(double *values, size_t len, double rate)
{
    for (size_t i = 0; i < len; i++) {
        values[i] = rate;
    }
}

void tally_tables_build(tally_tables *tables, const tally_config *config)
{
    phred_build(&tables->quality, config->min_phred);
    tables->uniform = config->uniform;
    tables->band    = config->band;
}

/* Points each array of rates at its own len values of storage. */
static void bind_arrays(phmm_rates *rates, const double *storage, size_t len)
{
    rates->modification   = storage;
    rates->open_insertion = storage + len;
    rates->open_deletion  = storage + 2 * len;
    rates->ends           = storage + 3 * len;
}

/* Writes the uniform rates into storage, in the order bind_arrays lays the arrays out. */
static void fill_arrays(double *storage, size_t len, const phmm_uniform_rates *uniform)
{
    fill_uniform(storage, len, uniform->modification);
    fill_uniform(storage + len, len, uniform->open_insertion);
    fill_uniform(storage + 2 * len, len, uniform->open_deletion);
    fill_uniform(storage + 3 * len, len, UNIFORM_ENDS);
}

void tally_rates_fill(phmm_rates *rates, double *storage, size_t len,
                      const tally_tables *tables)
{
    phmm_rates_set_transitions(rates, &tables->uniform);
    bind_arrays(rates, storage, len);
    fill_arrays(storage, len, &tables->uniform);
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

/* Returns the accumulator that counts this record once it is tallied. A supplementary
 * record places a piece of a read, and is counted apart from the read itself. */
static accum_field_id counted_field(const cm_bam_record *read)
{
    return filter_is_supplementary(read) ? ACCUM_SUPPLEMENTARY_COUNTED
                                         : ACCUM_PRIMARY_COUNTED;
}

/* Returns the accumulator that counts this record once it is turned away. */
static accum_field_id rejected_field(const cm_bam_record *read)
{
    return filter_is_supplementary(read) ? ACCUM_SUPPLEMENTARY_REJECTED
                                         : ACCUM_PRIMARY_REJECTED;
}

phmm_status tally(const cm_bam_record *read, const cm_fasta_record *ref,
                  const phmm_rates *rates, const tally_tables *tables,
                  tally_scratch *scratch, accum *target, pairs *target_pairs)
{
    context ctx = {
        .read         = read,
        .ref          = ref,
        .rates        = rates,
        .tables       = tables,
        .target       = target,
        .target_pairs = target_pairs,
    };
    phmm_status status = marginalize(&ctx, scratch);

    /* A record the model gives no path is counted where a record a filter turned away is
     * counted, and adds no other value. No value has reached the target: the window is added
     * only on PHMM_OK, and the failure is seen before either count below. */
    if (status == PHMM_NO_PATH) {
        *accum_data(target, rejected_field(read)) += 1.0;
        return status;
    }

    if (status != PHMM_OK) {
        return status;
    }

    *accum_data(target, counted_field(read)) += 1.0;

    add_length(&ctx);

    return PHMM_OK;
}
