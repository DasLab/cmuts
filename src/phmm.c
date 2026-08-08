/* phmm.c -- a banded pair HMM over one read and the reference it was placed on.
 *
 * Three states, the usual ones: M pairs a read base with a reference base, I
 * takes a read base the reference has nothing for, and D passes over a
 * reference base the read has nothing for. Both indices are prefix lengths, so
 * a cell (i, j) is the first i placed read bases set against the first j bases
 * of the reference, and the state says what the last operation was. A cell in
 * M or I therefore reads reference base j - 1 and read base i - 1.
 *
 * Only a band of reference positions is carried for each read base, centered on
 * where the CIGAR put it. Everything below is thus a probability conditioned on
 * the alignment staying inside that band, which is the approximation the whole
 * thing rests on: local ambiguity is marginalized away, and a grossly misplaced
 * read is not rescued.
 *
 * Rows are divided by their own totals as the forward pass goes and the
 * backward pass divides by the same numbers again, which is what lets a
 * posterior be read off as a product of the two with no separate normalizer and
 * no logarithms. Those divisors multiply back up to the likelihood of the read,
 * so nothing is lost by never forming it.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "phmm.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "align.h"
#include "nuc.h"

/* Reference positions carried for each read base. */
#define BAND_WIDTH (2 * PHMM_BAND + 1)

enum { STATE_MATCH, STATE_INSERTION, STATE_DELETION, N_STATES };

/* A base drawn from nothing in particular, and a comparison that settles
 * nothing: an inserted base answers to no reference position, and a base
 * neither side has named says as little either way. */
#define UNINFORMATIVE 0.25

/* How far the reference a read may reach before it is left to the walk. A band
 * following a CIGAR follows its reference skips too, and a read spanning an
 * intron would otherwise size its window to the intron. */
#define MAX_SPAN (PHMM_MAX_QUERY + 2 * PHMM_BAND + 2)

/* Forward times backward has to come to one on every row. A departure past this
 * is an index gone wrong rather than a rounding, and the read is better handed
 * back than counted from a matrix that does not hold. */
#define NORMALIZATION_TOLERANCE 1e-6

struct phmm_scratch {
    hts_pos_t *centers;    /* one per placed read base, and one before them */
    double    *forward;    /* rows, each BAND_WIDTH cells of N_STATES */
    double    *scale;      /* what each forward row was divided by */
    double    *coverage;   /* the window handed back */
    double    *spanned;
    double    *mutations;
    size_t     rows;       /* rows the three buffers above are sized for */
    size_t     window;     /* positions the two after them are sized for */

    /* Only the row in hand and the one below it are ever wanted, and the band
     * is a fixed width, so this one never has to grow. */
    double     backward[2][BAND_WIDTH][N_STATES];
};

/* Everything one read is marginalized against, gathered so that the passes read
 * as the recursions they are. */
typedef struct {
    const phmm            *model;
    const phred           *quality;
    const cm_bam_record   *read;
    const cm_fasta_record *ref;
    phmm_scratch          *scratch;
    aln_span               span;
    size_t                 rows;    /* placed bases, and one row before them */
    hts_pos_t              origin;  /* reference position of window value 0 */
    size_t                 window;
} context;

/* ------------------------------------------------------------------------ */
/* The model                                                                 */
/* ------------------------------------------------------------------------ */

phmm_params phmm_defaults(void)
{
    return (phmm_params){
        .open_deletion    = 1e-3,
        .open_insertion   = 2e-4,
        .extend_deletion  = 0.35,
        .extend_insertion = 0.35,
        .modification     = 5e-3,
    };
}

phmm_weights phmm_default_weights(void)
{
    return (phmm_weights){
        .weight = {
            [PHMM_SUBSTITUTION] = 1.0,
            [PHMM_DELETION]     = 1.0,
            [PHMM_INSERTION]    = 1.0,
        },
    };
}

void phmm_build(phmm *model, const phmm_params *params,
                const phmm_weights *weights)
{
    model->params                 = *params;
    model->weights                = *weights;
    model->match_to_insertion     = params->open_insertion;
    model->match_to_deletion      = params->open_deletion;
    model->match_to_match         = 1.0 - params->open_insertion
                                        - params->open_deletion;
    model->insertion_to_insertion = params->extend_insertion;
    model->insertion_to_match     = 1.0 - params->extend_insertion;
    model->deletion_to_deletion   = params->extend_deletion;
    model->deletion_to_match      = 1.0 - params->extend_deletion;
}

/* Two ways for a read base to disagree with the reference and one to agree with
 * it: the template may have been modified, and the read of it may be wrong.
 * Both are unlikely and neither is impossible, and it is the second alone that
 * a quality score speaks to. Agreement plus three disagreements comes to one,
 * which is what makes these a distribution over the four bases. */
static double agreement_chance(double modification, double error)
{
    return (1.0 - modification) * (1.0 - error) + modification * error / 3.0;
}

static double disagreement_chance(double modification, double error)
{
    return (1.0 - modification) * error / 3.0
         + modification * (1.0 - error / 3.0) / 3.0;
}

double phmm_modification(const phmm *model, bool agree, double error)
{
    double modification = model->params.modification;

    return agree
         ? modification * error / 3.0
             / agreement_chance(modification, error)
         : modification * (1.0 - error / 3.0) / 3.0
             / disagreement_chance(modification, error);
}

/* ------------------------------------------------------------------------ */
/* One read base against one reference base                                  */
/* ------------------------------------------------------------------------ */

/* A record storing no qualities weighs every base fully, there being nothing to
 * say otherwise, which leaves the modification rate alone to explain a
 * disagreement. */
static double error_at(const context *ctx, int32_t query)
{
    return ctx->read->qual
         ? phred_error(ctx->quality, ctx->read->qual[query])
         : 0.0;
}

typedef struct {
    bool   comparable;  /* both are named bases, and both are there to read */
    bool   agree;
    double error;
} comparison;

/* Cell (i, j) reads read base i - 1 of the placed span against reference base
 * j - 1. A position past either end of the reference has nothing to be compared
 * with, which is a thing not known rather than a disagreement. */
static comparison compare_at(const context *ctx, size_t i, hts_pos_t j)
{
    comparison out   = { .comparable = false, .agree = false, .error = 0.0 };
    int32_t    query = ctx->span.begin + (int32_t)i - 1;
    nuc        theirs;
    nuc        ours;

    if (j < 1 || (size_t)j > ctx->ref->len)
        return out;

    theirs = nuc_from_char(ctx->ref->seq[j - 1]);
    ours   = nuc_from_read(ctx->read->seq, query);

    if (!nuc_is_base(theirs) || !nuc_is_base(ours))
        return out;

    out.comparable = true;
    out.agree      = theirs == ours;
    out.error      = error_at(ctx, query);

    return out;
}

static double match_emission(const context *ctx, size_t i, hts_pos_t j)
{
    comparison at           = compare_at(ctx, i, j);
    double     modification = ctx->model->params.modification;

    if (!at.comparable)
        return UNINFORMATIVE;

    return at.agree ? agreement_chance(modification, at.error)
                    : disagreement_chance(modification, at.error);
}

static double modification_at(const context *ctx, size_t i, hts_pos_t j)
{
    comparison at = compare_at(ctx, i, j);

    return at.comparable ? phmm_modification(ctx->model, at.agree, at.error)
                         : 0.0;
}

/* How far the base a row pairs is to be believed, whatever it turned out to be
 * set against: a base neither side could name was still read, and read as well
 * or as badly as its score says. */
static double confidence_at(const context *ctx, size_t i)
{
    return 1.0 - error_at(ctx, ctx->span.begin + (int32_t)i - 1);
}

/* ------------------------------------------------------------------------ */
/* The band                                                                  */
/* ------------------------------------------------------------------------ */

static double *row_of(const context *ctx, size_t i)
{
    return ctx->scratch->forward + i * BAND_WIDTH * N_STATES;
}

/* The reference prefix length a cell of a row stands for. Every row is the same
 * width and moves with the CIGAR, so two rows are a shift apart and a cell on
 * one is found on the other by adding it. */
static hts_pos_t position_of(const context *ctx, size_t i, int k)
{
    return ctx->scratch->centers[i] - PHMM_BAND + k;
}

static hts_pos_t shift_between(const context *ctx, size_t from, size_t to)
{
    return ctx->scratch->centers[to] - ctx->scratch->centers[from];
}

static bool within_band(hts_pos_t k)
{
    return k >= 0 && k < BAND_WIDTH;
}

/* ------------------------------------------------------------------------ */
/* Forward                                                                   */
/* ------------------------------------------------------------------------ */

/* Nothing has been read yet, so each cell of the first row is the alignment
 * poised to begin at that reference position, none of them preferred over the
 * others. None has begun with a deletion: a reference base passed over before
 * the read starts is a base the read says nothing about, and its position could
 * never be recovered. */
static void forward_first_row(const context *ctx)
{
    double *row = row_of(ctx, 0);

    for (int k = 0; k < BAND_WIDTH; k++) {
        row[k * N_STATES + STATE_MATCH]     = 1.0 / BAND_WIDTH;
        row[k * N_STATES + STATE_INSERTION] = 0.0;
        row[k * N_STATES + STATE_DELETION]  = 0.0;
    }
}

/* The two states that take a read base, and so read from the row above. */
static void forward_paired_states(const context *ctx, size_t i)
{
    const phmm   *model = ctx->model;
    double       *row   = row_of(ctx, i);
    const double *above = row_of(ctx, i - 1);
    hts_pos_t     shift = shift_between(ctx, i - 1, i);

    for (int k = 0; k < BAND_WIDTH; k++) {
        hts_pos_t diagonal = k - 1 + shift;  /* cell for j - 1, one row up */
        hts_pos_t straight = k + shift;      /* cell for j,     one row up */
        double    paired   = 0.0;
        double    inserted = 0.0;

        if (within_band(diagonal)) {
            const double *from = above + diagonal * N_STATES;

            paired = model->match_to_match     * from[STATE_MATCH]
                   + model->insertion_to_match * from[STATE_INSERTION]
                   + model->deletion_to_match  * from[STATE_DELETION];
            paired *= match_emission(ctx, i, position_of(ctx, i, k));
        }

        if (within_band(straight)) {
            const double *from = above + straight * N_STATES;

            inserted = model->match_to_insertion     * from[STATE_MATCH]
                     + model->insertion_to_insertion * from[STATE_INSERTION];
            inserted *= UNINFORMATIVE;
        }

        row[k * N_STATES + STATE_MATCH]     = paired;
        row[k * N_STATES + STATE_INSERTION] = inserted;
    }
}

/* A deletion runs along the reference without the read moving, so it is the one
 * state whose row depends on itself: a cell takes from the cell to its left,
 * which must already hold this row's answer and not the last one's. */
static void forward_deletions(const context *ctx, size_t i)
{
    const phmm *model = ctx->model;
    double     *row   = row_of(ctx, i);

    row[STATE_DELETION] = 0.0;

    for (int k = 1; k < BAND_WIDTH; k++) {
        const double *left = row + (k - 1) * N_STATES;

        row[k * N_STATES + STATE_DELETION] =
              model->match_to_deletion    * left[STATE_MATCH]
            + model->deletion_to_deletion * left[STATE_DELETION];
    }
}

static void clear_deletions(const context *ctx, size_t i)
{
    double *row = row_of(ctx, i);

    for (int k = 0; k < BAND_WIDTH; k++)
        row[k * N_STATES + STATE_DELETION] = 0.0;
}

/* A deletion may neither open a read nor close one: reference passed over
 * before the first read base or after the last is reference whose position
 * could never be recovered, and admitting it would leave the read carrying an
 * event nothing in it points at.
 *
 * Both passes ask this rather than each forbidding it their own way. A step one
 * of them allows and the other does not is a disagreement about which paths
 * exist, and posteriors read off the two would not sum to one. */
static bool deletions_live(const context *ctx, size_t i)
{
    return i > 0 && i + 1 < ctx->rows;
}

/* Every row is divided by its own total, so that the numbers stay near one
 * however long the read: unscaled, a forward pass underflows a double within a
 * few hundred bases. */
static bool rescale(const context *ctx, size_t i)
{
    double *row   = row_of(ctx, i);
    double  total = 0.0;

    for (int c = 0; c < BAND_WIDTH * N_STATES; c++)
        total += row[c];

    if (!(total > 0.0) || !isfinite(total))
        return false;

    for (int c = 0; c < BAND_WIDTH * N_STATES; c++)
        row[c] /= total;

    ctx->scratch->scale[i] = total;

    return true;
}

static bool forward(const context *ctx)
{
    forward_first_row(ctx);

    if (!rescale(ctx, 0))
        return false;

    for (size_t i = 1; i < ctx->rows; i++) {
        forward_paired_states(ctx, i);

        if (deletions_live(ctx, i))
            forward_deletions(ctx, i);
        else
            clear_deletions(ctx, i);

        if (!rescale(ctx, i))
            return false;
    }

    return true;
}

/* ------------------------------------------------------------------------ */
/* Backward                                                                  */
/* ------------------------------------------------------------------------ */

/* The alignment ends on the last row, having paired or inserted its final base
 * but not having passed over a reference base to do it. */
static void backward_last_row(const context *ctx)
{
    double (*row)[N_STATES] = ctx->scratch->backward[(ctx->rows - 1) & 1];

    for (int k = 0; k < BAND_WIDTH; k++) {
        row[k][STATE_MATCH]     = 1.0;
        row[k][STATE_INSERTION] = 1.0;
        row[k][STATE_DELETION]  = 0.0;
    }
}

/* Everything reached from a cell, which for the two states that stay on this
 * row means the cell to its right. Only what crosses to the row below carries
 * that row's divisor, the deletion chain along this one having been divided
 * already as it was written. */
static void backward_row(const context *ctx, size_t i)
{
    const phmm *model = ctx->model;
    hts_pos_t   shift = shift_between(ctx, i, i + 1);
    double      below_scale = ctx->scratch->scale[i + 1];
    bool        live  = deletions_live(ctx, i);

    double       (*row)[N_STATES]   = ctx->scratch->backward[i & 1];
    const double (*below)[N_STATES] = ctx->scratch->backward[(i + 1) & 1];

    for (int k = BAND_WIDTH; k-- > 0; ) {
        hts_pos_t j        = position_of(ctx, i, k);
        hts_pos_t diagonal = k + 1 - shift;  /* cell for j + 1, one row down */
        hts_pos_t straight = k - shift;      /* cell for j,     one row down */
        double    paired   = 0.0;
        double    inserted = 0.0;
        double    deleted  = live && k + 1 < BAND_WIDTH
                           ? row[k + 1][STATE_DELETION]
                           : 0.0;

        if (within_band(diagonal))
            paired = match_emission(ctx, i + 1, j + 1)
                   * below[diagonal][STATE_MATCH] / below_scale;

        if (within_band(straight))
            inserted = UNINFORMATIVE
                     * below[straight][STATE_INSERTION] / below_scale;

        row[k][STATE_MATCH] = model->match_to_match     * paired
                            + model->match_to_insertion * inserted
                            + model->match_to_deletion  * deleted;

        row[k][STATE_INSERTION] = model->insertion_to_match     * paired
                                + model->insertion_to_insertion * inserted;

        row[k][STATE_DELETION] = live
                               ? model->deletion_to_match    * paired
                               + model->deletion_to_deletion * deleted
                               : 0.0;
    }
}

/* ------------------------------------------------------------------------ */
/* What a row says about the reference                                       */
/* ------------------------------------------------------------------------ */

static void add_at(const context *ctx, double *field, hts_pos_t position,
                   double value)
{
    hts_pos_t offset = position - ctx->origin;

    if (position < 0 || (size_t)position >= ctx->ref->len)
        return;

    if (offset >= 0 && (size_t)offset < ctx->window)
        field[offset] += value;
}

/* The chance a deletion ends here rather than at any of the other places it
 * could have: the run in hand, the step back out of it, and the pairing that
 * follows. That pairing lies on the row below, so the divisor which carried
 * this row into it has to be undone.
 *
 * Counted where the run ends rather than where it begins because reverse
 * transcription reads the template from its 3' end, so the last base a deletion
 * passes over is the first the enzyme met, and the likeliest to carry what
 * stopped it. The two are the marginals of one joint over runs -- summing
 * either over the runs it can belong to gives the same expected number of
 * events -- so the choice moves where a deletion is counted and not how much of
 * it there is. */
static double closed_deletion(const context *ctx, size_t i, int k)
{
    const double *forward_row = row_of(ctx, i);
    const double (*below)[N_STATES] = ctx->scratch->backward[(i + 1) & 1];
    hts_pos_t     j = position_of(ctx, i, k);
    hts_pos_t     diagonal;

    if (i + 1 >= ctx->rows)
        return 0.0;

    diagonal = k + 1 - shift_between(ctx, i, i + 1);

    if (!within_band(diagonal))
        return 0.0;

    return forward_row[k * N_STATES + STATE_DELETION]
         * ctx->model->deletion_to_match
         * match_emission(ctx, i + 1, j + 1)
         * below[diagonal][STATE_MATCH]
         / ctx->scratch->scale[i + 1];
}

/* The same for an insertion, whose pairing lies on the row above: the divisor
 * that carried that row into this one has to be undone. */
static double opened_insertion(const context *ctx, size_t i, int k)
{
    const double (*back)[N_STATES] = ctx->scratch->backward[i & 1];
    hts_pos_t      above;

    if (i == 0)
        return 0.0;

    above = k + shift_between(ctx, i - 1, i);

    if (!within_band(above))
        return 0.0;

    return row_of(ctx, i - 1)[above * N_STATES + STATE_MATCH]
         * ctx->model->match_to_insertion
         * UNINFORMATIVE
         * back[k][STATE_INSERTION]
         / ctx->scratch->scale[i];
}

/* A pairing spans the position it pairs with and covers it as far as the base
 * is to be believed, and carries whatever of its posterior belongs to the
 * template having really differed. A deletion spans every position it passes
 * over, the read having reached them all, but covers none of them, there being
 * no base read there to believe; it is a modification only where it opened,
 * since one adduct stops one reverse transcriptase once, whatever length of
 * reference it then skipped. An insertion is counted the same way, at the
 * position it sits before, and neither spans nor covers anything.
 *
 * Every event laid down is then worth what its kind is worth, the three saying
 * different amounts about a modification whatever the posterior says about the
 * events themselves. The mutations alone are weighed: what a position was
 * reached by does not bear on its having been reached.
 *
 * The first row is the alignment poised to begin and not yet begun. Its
 * posterior says where the read starts rather than what any base of it was set
 * against, and there is no base of it to ask about, so it contributes to
 * nothing. */
static void accumulate_row(const context *ctx, size_t i)
{
    phmm_scratch *scratch     = ctx->scratch;
    const phmm   *model       = ctx->model;
    const double *forward_row = row_of(ctx, i);
    const double (*back)[N_STATES] = scratch->backward[i & 1];
    double        confidence;

    if (i == 0)
        return;

    confidence = confidence_at(ctx, i);

    for (int k = 0; k < BAND_WIDTH; k++) {
        hts_pos_t j      = position_of(ctx, i, k);
        double    paired = forward_row[k * N_STATES + STATE_MATCH]
                         * back[k][STATE_MATCH];
        double    passed = forward_row[k * N_STATES + STATE_DELETION]
                         * back[k][STATE_DELETION];

        add_at(ctx, scratch->coverage, j - 1, paired * confidence);
        add_at(ctx, scratch->spanned, j - 1, paired + passed);
        add_at(ctx, scratch->mutations, j - 1,
               phmm_weigh(model, PHMM_SUBSTITUTION,
                          paired * modification_at(ctx, i, j)));
        add_at(ctx, scratch->mutations, j - 1,
               phmm_weigh(model, PHMM_DELETION, closed_deletion(ctx, i, k)));
        add_at(ctx, scratch->mutations, j,
               phmm_weigh(model, PHMM_INSERTION, opened_insertion(ctx, i, k)));
    }
}

/* Forward and backward have to describe the same set of paths, or what they
 * multiply out to are not posteriors at all. Every row coming to one is exactly
 * that statement, and it is cheap enough to make of one of them. */
static bool normalized(const context *ctx)
{
    const double *forward_row = row_of(ctx, 0);
    const double (*back)[N_STATES] = ctx->scratch->backward[0];
    double        total = 0.0;

    for (int k = 0; k < BAND_WIDTH; k++)
        for (int s = 0; s < N_STATES; s++)
            total += forward_row[k * N_STATES + s] * back[k][s];

    return fabs(total - 1.0) < NORMALIZATION_TOLERANCE;
}

static bool backward(const context *ctx)
{
    backward_last_row(ctx);
    accumulate_row(ctx, ctx->rows - 1);

    for (size_t i = ctx->rows - 1; i-- > 0; ) {
        backward_row(ctx, i);
        accumulate_row(ctx, i);
    }

    return normalized(ctx);
}

/* ------------------------------------------------------------------------ */
/* Scratch                                                                   */
/* ------------------------------------------------------------------------ */

phmm_scratch *phmm_scratch_create(void)
{
    return calloc(1, sizeof(phmm_scratch));
}

void phmm_scratch_destroy(phmm_scratch *scratch)
{
    if (!scratch)
        return;

    free(scratch->centers);
    free(scratch->forward);
    free(scratch->scale);
    free(scratch->coverage);
    free(scratch->spanned);
    free(scratch->mutations);
    free(scratch);
}

/* Each buffer that grew is kept whether or not its neighbours did, so that a
 * scratch failing to grow is still the scratch it was, just larger in places
 * than it needs to be. */
static int grow_rows(phmm_scratch *scratch, size_t rows)
{
    hts_pos_t *centers;
    double    *forward;
    double    *scale;

    if (rows <= scratch->rows)
        return 0;

    centers = realloc(scratch->centers, rows * sizeof *centers);
    forward = realloc(scratch->forward,
                      rows * BAND_WIDTH * N_STATES * sizeof *forward);
    scale   = realloc(scratch->scale, rows * sizeof *scale);

    if (centers)
        scratch->centers = centers;
    if (forward)
        scratch->forward = forward;
    if (scale)
        scratch->scale = scale;

    if (!centers || !forward || !scale)
        return -1;

    scratch->rows = rows;

    return 0;
}

static int grow_window(phmm_scratch *scratch, size_t window)
{
    double *coverage;
    double *spanned;
    double *mutations;

    if (window <= scratch->window)
        return 0;

    coverage  = realloc(scratch->coverage, window * sizeof *coverage);
    spanned   = realloc(scratch->spanned, window * sizeof *spanned);
    mutations = realloc(scratch->mutations, window * sizeof *mutations);

    if (coverage)
        scratch->coverage = coverage;
    if (spanned)
        scratch->spanned = spanned;
    if (mutations)
        scratch->mutations = mutations;

    if (!coverage || !spanned || !mutations)
        return -1;

    scratch->window = window;

    return 0;
}

/* ------------------------------------------------------------------------ */
/* One read                                                                  */
/* ------------------------------------------------------------------------ */

/* Every position a cell of any row can name, from the one before the first
 * band to the last the final band reaches. Reference prefix lengths only ever
 * grow, so the two ends are the first and last centers. */
static void size_window(context *ctx)
{
    const hts_pos_t *centers = ctx->scratch->centers;

    ctx->origin = centers[0] - PHMM_BAND - 1;
    ctx->window = (size_t)(centers[ctx->rows - 1] - centers[0])
                + 2 * PHMM_BAND + 2;
}

static bool prepare(context *ctx)
{
    const cm_bam_record *read = ctx->read;

    if (!read->seq || read->l_qseq <= 0 || read->l_qseq > PHMM_MAX_QUERY)
        return false;

    if (grow_rows(ctx->scratch, (size_t)read->l_qseq + 1) < 0)
        return false;

    ctx->span = aln_centers(read, ctx->scratch->centers);

    if (ctx->span.end <= ctx->span.begin)
        return false;

    ctx->rows = (size_t)(ctx->span.end - ctx->span.begin) + 1;

    size_window(ctx);

    if (ctx->window > MAX_SPAN || grow_window(ctx->scratch, ctx->window) < 0)
        return false;

    memset(ctx->scratch->coverage, 0,
           ctx->window * sizeof *ctx->scratch->coverage);
    memset(ctx->scratch->spanned, 0,
           ctx->window * sizeof *ctx->scratch->spanned);
    memset(ctx->scratch->mutations, 0,
           ctx->window * sizeof *ctx->scratch->mutations);

    return true;
}

bool phmm_run(const phmm *model, const phred *quality,
              const cm_bam_record *read, const cm_fasta_record *ref,
              phmm_scratch *scratch, phmm_window *out)
{
    context ctx = {
        .model   = model,
        .quality = quality,
        .read    = read,
        .ref     = ref,
        .scratch = scratch,
    };

    if (!prepare(&ctx) || !forward(&ctx) || !backward(&ctx))
        return false;

    out->origin    = ctx.origin;
    out->len       = ctx.window;
    out->coverage  = scratch->coverage;
    out->spanned   = scratch->spanned;
    out->mutations = scratch->mutations;

    return true;
}
