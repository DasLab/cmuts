/* phmm.c -- a banded pair HMM over one read and the reference it was placed on.
 *
 * Three states, the usual ones: M pairs a read base with a reference base, I
 * takes a read base the reference has nothing for, and D passes over a
 * reference base the read has nothing for. Both indices are prefix lengths, so
 * a cell (i, j) is the first i placed read bases set against the first j bases
 * of the reference, and the state says what the last operation was. A cell in
 * M or I therefore reads reference base j - 1 and read base i - 1.
 *
 * Only a band of reference positions is carried for each read base: what the
 * CIGAR path crosses on that base's row, and a half-width either side of it.
 * Everything below is thus a probability conditioned on the alignment staying
 * inside that band, which is the approximation the whole thing rests on: local
 * ambiguity is marginalized away, and a grossly misplaced read is not rescued.
 * A band of nothing leaves the CIGAR path alone and marginalizes over it.
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

enum { STATE_MATCH, STATE_INSERTION, STATE_DELETION, N_STATES };

/* One cell of a band row, named so that a row of them can be returned. */
typedef double band_cell[N_STATES];

/* A base drawn from nothing in particular, and a comparison that settles
 * nothing: an inserted base answers to no reference position, and a base
 * neither side has named says as little either way. */
#define UNINFORMATIVE 0.25

/* Forward times backward has to come to one on every row. A departure past this
 * is an index gone wrong rather than a rounding, and the read is better handed
 * back than counted from a matrix that does not hold. */
#define NORMALIZATION_TOLERANCE 1e-6

/* Rows are stored at the widest row's stride rather than each at its own, so
 * that a row is found by multiplication and no table of offsets is needed. It
 * is the loops that are ragged, and they are what the cost follows. */
struct phmm_scratch {
    aln_place *places;     /* one per placed read base, and one before them */
    double    *scale;      /* what each forward row was divided by */
    double    *forward;    /* rows, each of widest cells of N_STATES */
    /* Only the row in hand and the one below it are ever wanted, so two are
     * kept however many rows the read has. */
    band_cell *backward;
    double    *coverage;   /* the window handed back */
    double    *spanned;
    double    *mutations;
    size_t     rows;         /* rows places and scale are sized for */
    size_t     matrix_rows;  /* rows the forward matrix is sized for */
    size_t     widest;       /* cells a row of it holds */
    size_t     window;       /* positions the last three are sized for */
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
    const int             *half;    /* how far either side of the CIGAR each row
                                       may look; the caller's, one per row */
    hts_pos_t              widest;  /* cells the widest of them holds */
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
    return ctx->scratch->forward + i * (size_t)ctx->widest * N_STATES;
}

/* Which of the two a row lands in is its parity. */
static band_cell *backward_row_of(const context *ctx, size_t i)
{
    return ctx->scratch->backward + (i & 1) * (size_t)ctx->widest;
}

/* The reference a row's own base accounts for beyond the one it pairs with,
 * which is however much a deletion after it passes over. */
static hts_pos_t skip_at(const context *ctx, size_t i)
{
    return ctx->scratch->places[i].last - ctx->scratch->places[i].first;
}

/* A row covers everything the CIGAR path crosses on it and its half-width
 * either side of that, rather than a half-width either side of a point. The
 * distinction is the deletion: the path runs along a row to cross one, so a row
 * with a deletion after it is a stretch that a band of nothing still holds
 * whole. A read is never denied the alignment it arrived with, and the band
 * means room to depart from it and nothing else. */
static hts_pos_t width_at(const context *ctx, size_t i)
{
    return skip_at(ctx, i) + 2 * (hts_pos_t)ctx->half[i] + 1;
}

/* The reference prefix length a row's first cell stands for. */
static hts_pos_t origin_of(const context *ctx, size_t i)
{
    return ctx->scratch->places[i].first - ctx->half[i];
}

/* The reference prefix length a cell of a row stands for. */
static hts_pos_t position_of(const context *ctx, size_t i, hts_pos_t k)
{
    return origin_of(ctx, i) + k;
}

/* What to add to a cell index to name the same position on another row. Rows
 * move with the CIGAR and may widen or narrow as they go, so what separates two
 * of them is the distance between their first cells. */
static hts_pos_t shift_between(const context *ctx, size_t from, size_t to)
{
    return origin_of(ctx, to) - origin_of(ctx, from);
}

/* Rows no longer share a width, so a cell index means nothing apart from the
 * row it indexes. Asking it of the wrong row is how the two passes would come
 * to disagree about which paths exist. */
static bool within_row(const context *ctx, size_t i, hts_pos_t k)
{
    return k >= 0 && k < width_at(ctx, i);
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
    double *row   = row_of(ctx, 0);
    hts_pos_t width = width_at(ctx, 0);

    for (int k = 0; k < width; k++) {
        row[k * N_STATES + STATE_MATCH]     = 1.0 / width;
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

    for (hts_pos_t k = 0; k < width_at(ctx, i); k++) {
        hts_pos_t diagonal = k - 1 + shift;  /* cell for j - 1, one row up */
        hts_pos_t straight = k + shift;      /* cell for j,     one row up */
        double    paired   = 0.0;
        double    inserted = 0.0;

        if (within_row(ctx, i - 1, diagonal)) {
            const double *from = above + diagonal * N_STATES;

            paired = model->match_to_match     * from[STATE_MATCH]
                   + model->insertion_to_match * from[STATE_INSERTION]
                   + model->deletion_to_match  * from[STATE_DELETION];
            paired *= match_emission(ctx, i, position_of(ctx, i, k));
        }

        if (within_row(ctx, i - 1, straight)) {
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

    for (hts_pos_t k = 1; k < width_at(ctx, i); k++) {
        const double *left = row + (k - 1) * N_STATES;

        row[k * N_STATES + STATE_DELETION] =
              model->match_to_deletion    * left[STATE_MATCH]
            + model->deletion_to_deletion * left[STATE_DELETION];
    }
}

static void clear_deletions(const context *ctx, size_t i)
{
    double *row = row_of(ctx, i);

    for (hts_pos_t k = 0; k < width_at(ctx, i); k++)
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
    hts_pos_t cells = width_at(ctx, i) * N_STATES;
    double  total = 0.0;

    for (hts_pos_t c = 0; c < cells; c++)
        total += row[c];

    if (!(total > 0.0) || !isfinite(total))
        return false;

    for (hts_pos_t c = 0; c < cells; c++)
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
    band_cell *row = backward_row_of(ctx, ctx->rows - 1);

    for (hts_pos_t k = 0; k < width_at(ctx, ctx->rows - 1); k++) {
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

    band_cell       *row   = backward_row_of(ctx, i);
    const band_cell *below = backward_row_of(ctx, i + 1);

    for (hts_pos_t k = width_at(ctx, i); k-- > 0; ) {
        hts_pos_t j        = position_of(ctx, i, k);
        hts_pos_t diagonal = k + 1 - shift;  /* cell for j + 1, one row down */
        hts_pos_t straight = k - shift;      /* cell for j,     one row down */
        double    paired   = 0.0;
        double    inserted = 0.0;
        double    deleted  = live && k + 1 < width_at(ctx, i)
                           ? row[k + 1][STATE_DELETION]
                           : 0.0;

        if (within_row(ctx, i + 1, diagonal))
            paired = match_emission(ctx, i + 1, j + 1)
                   * below[diagonal][STATE_MATCH] / below_scale;

        if (within_row(ctx, i + 1, straight))
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
static double closed_deletion(const context *ctx, size_t i, hts_pos_t k)
{
    const double    *forward_row = row_of(ctx, i);
    const band_cell *below       = backward_row_of(ctx, i + 1);
    hts_pos_t        j           = position_of(ctx, i, k);
    hts_pos_t        diagonal;

    if (i + 1 >= ctx->rows)
        return 0.0;

    diagonal = k + 1 - shift_between(ctx, i, i + 1);

    if (!within_row(ctx, i + 1, diagonal))
        return 0.0;

    return forward_row[k * N_STATES + STATE_DELETION]
         * ctx->model->deletion_to_match
         * match_emission(ctx, i + 1, j + 1)
         * below[diagonal][STATE_MATCH]
         / ctx->scratch->scale[i + 1];
}

/* The same for an insertion, whose pairing lies on the row above: the divisor
 * that carried that row into this one has to be undone. */
static double opened_insertion(const context *ctx, size_t i, hts_pos_t k)
{
    const band_cell *back = backward_row_of(ctx, i);
    hts_pos_t        above;

    if (i == 0)
        return 0.0;

    above = k + shift_between(ctx, i - 1, i);

    if (!within_row(ctx, i - 1, above))
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
    phmm_scratch    *scratch     = ctx->scratch;
    const phmm      *model       = ctx->model;
    const double    *forward_row = row_of(ctx, i);
    const band_cell *back        = backward_row_of(ctx, i);
    double           confidence;

    if (i == 0)
        return;

    confidence = confidence_at(ctx, i);

    for (hts_pos_t k = 0; k < width_at(ctx, i); k++) {
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
    const double    *forward_row = row_of(ctx, 0);
    const band_cell *back        = backward_row_of(ctx, 0);
    double           total       = 0.0;

    for (hts_pos_t k = 0; k < width_at(ctx, 0); k++)
        for (hts_pos_t s = 0; s < N_STATES; s++)
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

    free(scratch->places);
    free(scratch->forward);
    free(scratch->scale);
    free(scratch->backward);
    free(scratch->coverage);
    free(scratch->spanned);
    free(scratch->mutations);
    free(scratch);
}

/* Each buffer that grew is kept whether or not its neighbours did, so that a
 * scratch failing to grow is still the scratch it was, just larger in places
 * than it needs to be. Nothing here ever shrinks, a scratch being reused across
 * every read one worker sees.
 *
 * The rows come first and alone, because how many there are is read off the
 * CIGAR, and the places must be somewhere to be written before they can say
 * how wide the widest row is. */
static int grow_rows(phmm_scratch *scratch, size_t rows)
{
    aln_place *places;
    double    *scale;

    if (rows <= scratch->rows)
        return 0;

    places = realloc(scratch->places, rows * sizeof *places);
    scale  = realloc(scratch->scale, rows * sizeof *scale);

    if (places)
        scratch->places = places;
    if (scale)
        scratch->scale = scale;

    if (!places || !scale)
        return -1;

    scratch->rows = rows;

    return 0;
}

/* Sized by the rows and the widest row together, so a band wider than the last
 * read's regrows them however few rows this one needs. */
static int grow_band(phmm_scratch *scratch, size_t rows, size_t widest)
{
    double    *forward;
    band_cell *backward;

    if (rows <= scratch->matrix_rows && widest <= scratch->widest)
        return 0;

    rows   = rows   > scratch->matrix_rows ? rows   : scratch->matrix_rows;
    widest = widest > scratch->widest      ? widest : scratch->widest;

    forward  = realloc(scratch->forward,
                       rows * widest * N_STATES * sizeof *forward);
    backward = realloc(scratch->backward, 2 * widest * sizeof *backward);

    if (forward)
        scratch->forward = forward;
    if (backward)
        scratch->backward = backward;

    if (!forward || !backward)
        return -1;

    scratch->matrix_rows = rows;
    scratch->widest      = widest;

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

/* The widest row, which is the stride every row is stored at and the only width
 * the buffers have to answer for. One at least, so that a band asked for in bad
 * faith leaves a matrix that can be addressed rather than one that cannot; the
 * rows of it are then empty and the pass dies on the first of them. */
static hts_pos_t widest_row(const context *ctx)
{
    hts_pos_t widest = 1;

    for (size_t i = 0; i < ctx->rows; i++)
        if (width_at(ctx, i) > widest)
            widest = width_at(ctx, i);

    return widest;
}

/* Every position a cell of any row can name, from the one before the earliest
 * a band reaches to the last one does. A row's reach is its center give or take
 * its own half-width, and a row further along may be narrower than one behind
 * it, so neither end belongs to a particular row and both are looked for. */
static void size_window(context *ctx)
{
    hts_pos_t first = origin_of(ctx, 0);
    hts_pos_t last  = first + width_at(ctx, 0);

    for (size_t i = 1; i < ctx->rows; i++) {
        hts_pos_t lo = origin_of(ctx, i);
        hts_pos_t hi = lo + width_at(ctx, i);

        if (lo < first)
            first = lo;
        if (hi > last)
            last = hi;
    }

    ctx->origin = first - 1;
    ctx->window = (size_t)(last - first) + 1;
}

static bool prepare(context *ctx)
{
    const cm_bam_record *read = ctx->read;

    if (grow_rows(ctx->scratch, (size_t)read->l_qseq + 1) < 0)
        return false;

    ctx->span   = aln_places(read, ctx->scratch->places);
    ctx->rows   = (size_t)(ctx->span.end - ctx->span.begin) + 1;
    ctx->widest = widest_row(ctx);

    if (grow_band(ctx->scratch, ctx->rows, (size_t)ctx->widest) < 0)
        return false;

    size_window(ctx);

    if (grow_window(ctx->scratch, ctx->window) < 0)
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
              const int *half, phmm_scratch *scratch, phmm_window *out)
{
    context ctx = {
        .model   = model,
        .quality = quality,
        .read    = read,
        .ref     = ref,
        .scratch = scratch,
        .half    = half,
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
