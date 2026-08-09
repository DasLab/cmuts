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

/* Everything a cell's comparison is afterwards asked for.
 *
 * The forward pass wants what the pairing is worth, the backward pass wants the
 * same of a row it has yet to reach, the deletion closing there wants it again,
 * and the accumulation wants how much of the pairing was a modification. All
 * four come of one comparison, so it is made once and held rather than made
 * again wherever it is wanted. Kept for the whole matrix rather than a row at
 * a time: the backward pass reads a row the forward pass wrote and left
 * behind. */
typedef struct {
    double emission;
    double modification;
} cell_terms;

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
    double    *scale;      /* what each forward row is scaled by, which is the
                              reciprocal of what it came to */
    band_cell *forward;    /* rows, each of widest cells */
    cell_terms *terms;     /* one per cell of it */
    /* Only the row in hand and the one below it are ever wanted, so two are
     * kept however many rows the read has. */
    band_cell *backward;
    /* What a pairing on the row below is worth, one per cell of the row in
     * hand. The backward pass forms it and the deletion closing there is a
     * product with it, so it is handed along rather than made twice. */
    double    *pairings;
    double    *coverage;   /* the window handed back */
    double    *spanned;
    double    *mutations;
    size_t     rows;         /* rows places and scale are sized for */
    size_t     matrix_rows;  /* rows the forward matrix is sized for */
    size_t     widest;       /* cells a row of it holds */
    size_t     window;       /* positions the last three are sized for */
};

/* A stretch of the reference, named by where it begins and how far it runs. */
typedef struct {
    hts_pos_t origin;
    size_t    len;
} extent;

/* Everything one read is marginalized against, gathered so that the passes read
 * as the recursions they are.
 *
 * What a caller hands over comes first and is fixed for the read. What follows
 * is worked out from it, and stands for nothing until prepare has run: a pass
 * reading it before then reads the band of whichever read came before. prepare
 * is the only thing that writes any of it, which is why every other function
 * here takes a context it may not modify. */
typedef struct {
    const phmm            *model;
    const phred           *quality;
    const cm_bam_record   *read;
    const cm_fasta_record *ref;
    phmm_scratch          *scratch;
    const int             *half;    /* how far either side of the CIGAR each row
                                       may look; the caller's, one per row */

    aln_span               span;    /* the stretch of the read that is placed */
    size_t                 rows;    /* placed bases, and one row before them */
    hts_pos_t              widest;  /* cells the widest row holds */
    extent                 window;  /* every position a row lays anything at */
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

/* Having read what was read, the chance the template really differed from the
 * reference here.
 *
 * A base agreeing with the reference may have been modified and then misread
 * back into agreement, and one disagreeing may be an unmodified base misread.
 * The first is negligible and the second is not, which is what makes a poorly
 * read disagreement worth less than a clean one. */
static double phmm_modification(const phmm *model, bool agree, double error)
{
    double modification = model->params.modification;

    return agree
         ? modification * error / 3.0
             / agreement_chance(modification, error)
         : modification * (1.0 - error / 3.0) / 3.0
             / disagreement_chance(modification, error);
}

/* ------------------------------------------------------------------------ */
/* One read base against the reference                                       */
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

/* What a pairing that agreed, or one that did not, is worth, and how much of it
 * belongs to the template having really differed. */
static cell_terms terms_from(const context *ctx, bool agree, double error)
{
    double m = ctx->model->params.modification;

    return (cell_terms){
        .emission     = agree ? agreement_chance(m, error)
                              : disagreement_chance(m, error),
        .modification = phmm_modification(ctx->model, agree, error),
    };
}

/* Every value a cell of one row can take.
 *
 * A row pairs one read base, which is read once and scored once, so the only
 * thing a cell of it settles is which of three the reference offers: the same
 * base, a different one, or nothing that can be compared at all. The three are
 * worked out when the row is entered and a cell picks among them, which is
 * what keeps the divisions a modification costs off the cell.
 *
 * A read that names no base here agrees with nothing and differs from nothing,
 * so all three are the same and a cell need not ask about the read again. */
typedef struct {
    cell_terms agree;
    cell_terms differ;
    cell_terms neither;
    nuc        ours;
} row_terms;

static row_terms row_terms_of(const context *ctx, size_t i)
{
    int32_t    query   = ctx->span.begin + (int32_t)i - 1;
    double     error   = error_at(ctx, query);
    nuc        ours    = nuc_from_read(ctx->read->seq, query);
    cell_terms neither = { .emission = UNINFORMATIVE, .modification = 0.0 };
    bool       named   = nuc_is_base(ours);

    return (row_terms){
        .agree   = named ? terms_from(ctx, true, error) : neither,
        .differ  = named ? terms_from(ctx, false, error) : neither,
        .neither = neither,
        .ours    = ours,
    };
}

/* Cell (i, j) reads the row's base against reference base j - 1. A position
 * past either end of the reference has nothing to be compared with, which is a
 * thing not known rather than a disagreement. */
static cell_terms terms_at(const context *ctx, const row_terms *row,
                           hts_pos_t j)
{
    nuc theirs;

    if (j < 1 || (size_t)j > ctx->ref->len)
        return row->neither;

    theirs = nuc_from_char(ctx->ref->seq[j - 1]);

    if (!nuc_is_base(theirs))
        return row->neither;

    return theirs == row->ours ? row->agree : row->differ;
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

/* A forward row as the pass left it, which is not what the row stands for: see
 * scaled_row below before reading one. */
static band_cell *row_of(const context *ctx, size_t i)
{
    return ctx->scratch->forward + i * (size_t)ctx->widest;
}

/* A forward row together with what restores it.
 *
 * A row is left scaled by every divisor but its own, that last one having been
 * folded into the transitions leaving the row rather than walked over the row
 * again. So what is stored is the row's own total times what the row stands
 * for, and the two are a probability only together.
 *
 * Undoing it belongs here and nowhere else. Reading a forward row any other way
 * is the one mistake this file cannot catch: the numbers look like posteriors
 * and are wrong by a factor that changes from row to row, so nothing overflows,
 * nothing fails to sum, and the counts come out quietly wrong.
 *
 * A row and its divisor are fetched together, once, because a pass that laid
 * anything down between fetching them would have to fetch them both again: the
 * window it writes is doubles and so is a row, and no compiler may assume the
 * two do not overlap. */
typedef struct {
    const band_cell *cell;
    double           scale;
} scaled_row;

static scaled_row scaled_row_of(const context *ctx, size_t i)
{
    return (scaled_row){
        .cell  = row_of(ctx, i),
        .scale = ctx->scratch->scale[i],
    };
}

/* What one cell of it stands for. */
static double forward_at(const scaled_row *row, hts_pos_t k, int state)
{
    return row->cell[k][state] * row->scale;
}

static cell_terms *terms_of(const context *ctx, size_t i)
{
    return ctx->scratch->terms + i * (size_t)ctx->widest;
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

/* Rows do not share a width, so a cell index means nothing apart from the row
 * it indexes, and it is that row's width the index has to be held against.
 * Holding it against another's is how the two passes would come to disagree
 * about which paths exist.
 *
 * The width is handed in rather than looked up because a row does not change
 * width while it is being walked, and fetching it again for every cell of every
 * neighbour is most of what the test costs. */
static bool within(hts_pos_t k, hts_pos_t width)
{
    return k >= 0 && k < width;
}

/* ------------------------------------------------------------------------ */
/* Forward                                                                   */
/* ------------------------------------------------------------------------ */

/* Nothing has been read yet, so each cell of the first row is the alignment
 * poised to begin at that reference position, none of them preferred over the
 * others. None has begun with a deletion: a reference base passed over before
 * the read starts is a base the read says nothing about, and its position could
 * never be recovered. */
static double forward_first_row(const context *ctx)
{
    band_cell *row   = row_of(ctx, 0);
    hts_pos_t  width = width_at(ctx, 0);
    double     share = 1.0 / (double)width;
    double     total = 0.0;

    for (hts_pos_t k = 0; k < width; k++) {
        row[k][STATE_MATCH]     = share;
        row[k][STATE_INSERTION] = 0.0;
        row[k][STATE_DELETION]  = 0.0;

        total += share;
    }

    return total;
}

/* What each state of the row above contributes to a cell of this one, with
 * everything constant along the row multiplied out once: what the row above is
 * scaled by, and the flat emission an inserted base carries. What is left is
 * the pairing's own emission, that being the one thing a cell decides for
 * itself.
 *
 * Folding the scaling in here is what spares the row a second walk. Dividing a
 * row by its total and stepping from it undivided leave the same numbers, but
 * the first costs a multiplication, a load and a store for every state of every
 * cell where this costs a handful of multiplications for the whole row.
 *
 * This is the one place that reads a forward row without going through
 * forward_at, and it is the same undoing by another route: a divisor that would
 * multiply every cell of the row above is carried instead by the five constants
 * every cell of it is met with. Everything else that reads a row goes through
 * forward_at, where what the undoing is for is written down. */
typedef struct {
    double match_to_match;
    double insertion_to_match;
    double deletion_to_match;
    double match_to_insertion;
    double insertion_to_insertion;
} descent;

static descent descent_into(const context *ctx, size_t i)
{
    const phmm *model = ctx->model;
    double      scale = ctx->scratch->scale[i - 1];

    return (descent){
        .match_to_match         = model->match_to_match * scale,
        .insertion_to_match     = model->insertion_to_match * scale,
        .deletion_to_match      = model->deletion_to_match * scale,
        .match_to_insertion     = model->match_to_insertion * scale
                                * UNINFORMATIVE,
        .insertion_to_insertion = model->insertion_to_insertion * scale
                                * UNINFORMATIVE,
    };
}

/* The two states that take a read base, and so read from the row above. */
static double paired_from(const descent *step, const double *above,
                          double emission)
{
    return (step->match_to_match     * above[STATE_MATCH]
          + step->insertion_to_match * above[STATE_INSERTION]
          + step->deletion_to_match  * above[STATE_DELETION]) * emission;
}

static double inserted_from(const descent *step, const double *above)
{
    return step->match_to_insertion     * above[STATE_MATCH]
         + step->insertion_to_insertion * above[STATE_INSERTION];
}

/* A deletion runs along the reference without the read moving, so it is the one
 * state whose row depends on itself: a cell takes from the cell to its left,
 * which must already hold this row's answer and not the last one's. */
static double deleted_from(const phmm *model, const double *left)
{
    return model->match_to_deletion    * left[STATE_MATCH]
         + model->deletion_to_deletion * left[STATE_DELETION];
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

/* One row of the forward pass, returning what it came to before it was scaled.
 *
 * Everything a cell wants is at hand by the time it is reached: its comparison
 * is made here, the two states that take a read base read the row above, and
 * the deletion reads the cell to the left, which this same loop wrote an
 * iteration ago. So it is one pass over the row rather than four, which at a
 * narrow band is most of what a row costs. The total is gathered as it goes,
 * and in the order the cells lie, so that the sum is the one the scaling
 * afterwards divides by. */
static double forward_row(const context *ctx, size_t i)
{
    const phmm      *model = ctx->model;
    descent          step  = descent_into(ctx, i);
    band_cell       *row   = row_of(ctx, i);
    const band_cell *above = row_of(ctx, i - 1);
    cell_terms      *terms = terms_of(ctx, i);
    row_terms        each  = row_terms_of(ctx, i);
    hts_pos_t        shift = shift_between(ctx, i - 1, i);
    hts_pos_t        width = width_at(ctx, i);
    hts_pos_t        above_width = width_at(ctx, i - 1);
    bool             live  = deletions_live(ctx, i);
    double           total = 0.0;

    for (hts_pos_t k = 0; k < width; k++) {
        hts_pos_t diagonal = k - 1 + shift;  /* a position back, one row up */
        hts_pos_t straight = k + shift;      /* this position, one row up */

        terms[k] = terms_at(ctx, &each, position_of(ctx, i, k));

        row[k][STATE_MATCH] = within(diagonal, above_width)
                            ? paired_from(&step, above[diagonal],
                                          terms[k].emission)
                            : 0.0;
        row[k][STATE_INSERTION] = within(straight, above_width)
                                ? inserted_from(&step, above[straight])
                                : 0.0;
        row[k][STATE_DELETION] = live && k > 0
                               ? deleted_from(model, row[k - 1])
                               : 0.0;

        total += row[k][STATE_MATCH];
        total += row[k][STATE_INSERTION];
        total += row[k][STATE_DELETION];
    }

    return total;
}

/* Every row is scaled by its own total, so that the numbers stay near one
 * however long the read: unscaled, a forward pass underflows a double within a
 * few hundred bases. The scaling itself is left to the descent out of the row,
 * so all that is kept here is what the descent will want. */
static bool record_total(const context *ctx, size_t i, double total)
{
    if (!(total > 0.0) || !isfinite(total))
        return false;

    ctx->scratch->scale[i] = 1.0 / total;

    return true;
}

/* The first row is laid out rather than stepped into, having no row above it. */
static bool forward(const context *ctx)
{
    if (!record_total(ctx, 0, forward_first_row(ctx)))
        return false;

    for (size_t i = 1; i < ctx->rows; i++)
        if (!record_total(ctx, i, forward_row(ctx, i)))
            return false;

    return true;
}

/* ------------------------------------------------------------------------ */
/* Backward                                                                  */
/* ------------------------------------------------------------------------ */

/* The alignment ends on the last row, having paired or inserted its final base
 * but not having passed over a reference base to do it. There is no row below
 * for a deletion to close on, so nothing is worth anything there. */
static void backward_last_row(const context *ctx)
{
    band_cell *row      = backward_row_of(ctx, ctx->rows - 1);
    double    *pairings = ctx->scratch->pairings;

    for (hts_pos_t k = 0; k < width_at(ctx, ctx->rows - 1); k++) {
        row[k][STATE_MATCH]     = 1.0;
        row[k][STATE_INSERTION] = 1.0;
        row[k][STATE_DELETION]  = 0.0;

        pairings[k] = 0.0;
    }
}

/* Everything reached from a cell, which for the two states that stay on this
 * row means the cell to its right. Only what crosses to the row below carries
 * that row's scaling, the deletion chain along this one having been scaled
 * already as it was written.
 *
 * What a pairing on the row below is worth is left behind for the accumulation,
 * the deletion closing there being a product with exactly this and nothing
 * else. */
static void backward_row(const context *ctx, size_t i)
{
    const phmm *model = ctx->model;
    hts_pos_t   shift = shift_between(ctx, i, i + 1);
    double      below_scale = ctx->scratch->scale[i + 1];
    bool        live  = deletions_live(ctx, i);

    band_cell        *row      = backward_row_of(ctx, i);
    const band_cell  *below    = backward_row_of(ctx, i + 1);
    const cell_terms *terms    = terms_of(ctx, i + 1);
    double           *pairings = ctx->scratch->pairings;
    hts_pos_t         width    = width_at(ctx, i);
    hts_pos_t         below_width = width_at(ctx, i + 1);

    for (hts_pos_t k = width; k-- > 0; ) {
        hts_pos_t diagonal = k + 1 - shift;  /* a position on, one row down */
        hts_pos_t straight = k - shift;      /* this position, one row down */
        double    paired   = 0.0;
        double    inserted = 0.0;
        double    deleted  = live && k + 1 < width
                           ? row[k + 1][STATE_DELETION]
                           : 0.0;

        if (within(diagonal, below_width))
            paired = terms[diagonal].emission
                   * below[diagonal][STATE_MATCH] * below_scale;

        if (within(straight, below_width))
            inserted = UNINFORMATIVE
                     * below[straight][STATE_INSERTION] * below_scale;

        pairings[k] = paired;

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

/* The three fields of the window, each advanced to where a row's first cell
 * enters it, so that a cell is addressed by its own index and no position is
 * held against the window's bounds.
 *
 * A cell at k stands for reference prefix length j: it pairs the base before
 * it and lays an insertion at j itself, so a row enters the window over
 * [j0 - 1, j0 + width - 1] and nowhere else. Every one of those is a position
 * the window holds, window_of running it from one before the earliest position
 * any row names to the furthest any row names. A row is therefore a stretch of
 * the window at a fixed offset, and that is an invariant of window_of rather
 * than of anything here: a window covering less would not admit this.
 *
 * A position past either end of the reference is not turned away. Some of what
 * the band names really does lie outside, a read placed near a boundary having
 * paths that leave it; which of them are worth keeping is the caller's to say,
 * and saying it twice would leave two places to look when the answer changed. */
typedef struct {
    double *coverage;
    double *spanned;
    double *mutations;
} landing;

static landing landing_of(const context *ctx, size_t i)
{
    phmm_scratch *scratch = ctx->scratch;
    size_t        at      = (size_t)(origin_of(ctx, i) - 1
                                   - ctx->window.origin);

    return (landing){
        .coverage  = scratch->coverage + at,
        .spanned   = scratch->spanned + at,
        .mutations = scratch->mutations + at,
    };
}

/* What one event of each kind, believed to whatever degree a cell believes it,
 * is worth to the position it is laid at.
 *
 * Each is the weight its kind carries and then everything the event's posterior
 * is a product with that a cell does not decide: for a deletion the step back
 * out of the run, and for an insertion the step into it and the flat emission
 * it carries. What is left for a cell is what a cell alone knows.
 *
 * An insertion also carries the divisor of the row it opens on, every posterior
 * over a step between two rows doing so. That is the recursion and not the way
 * a row is stored -- what a row is stored as is undone by forward_at, and the
 * two would otherwise be easy to read as one thing.
 *
 * Wanted only of a row that has one above it, an insertion opening on the first
 * row being an insertion nothing carried into. */
typedef struct {
    double substitution;
    double deletion;
    double insertion;
} weighing;

static weighing weighing_of(const context *ctx, size_t i)
{
    const phmm   *model  = ctx->model;
    const double *weight = model->weights.weight;

    return (weighing){
        .substitution = weight[PHMM_SUBSTITUTION],
        .deletion     = weight[PHMM_DELETION] * model->deletion_to_match,
        .insertion    = weight[PHMM_INSERTION] * model->match_to_insertion
                      * UNINFORMATIVE * ctx->scratch->scale[i],
    };
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
 * A deletion is laid where its run ends rather than where it begins because
 * reverse transcription reads the template from its 3' end, so the last base a
 * deletion passes over is the first the enzyme met, and the likeliest to carry
 * what stopped it. The two are the marginals of one joint over runs -- summing
 * either over the runs it can belong to gives the same expected number of
 * events -- so the choice moves where a deletion is counted and not how much of
 * it there is. What it ends against is the pairing on the row below, which the
 * backward pass left behind; what an insertion opens out of is the pairing on
 * the row above, which is still where the forward pass wrote it.
 *
 * The row's own scaling is undone on the two states that are asked about, that
 * being cheaper than undoing it on the whole row as it was written.
 *
 * The first row is the alignment poised to begin and not yet begun. Its
 * posterior says where the read starts rather than what any base of it was set
 * against, and there is no base of it to ask about, so it contributes to
 * nothing. */
static void accumulate_row(const context *ctx, size_t i)
{
    phmm_scratch     *scratch  = ctx->scratch;
    scaled_row        front    = scaled_row_of(ctx, i);
    const band_cell  *back     = backward_row_of(ctx, i);
    const cell_terms *terms    = terms_of(ctx, i);
    const double     *pairings = scratch->pairings;
    hts_pos_t         width    = width_at(ctx, i);
    /* An insertion belongs at the position its cell sits before, which is where
     * the next cell lays what it pairs, so it is held here and laid on the
     * iteration after. That leaves every cell entering one slot of the mutations
     * and no other, where laying it ahead would have each cell read back a slot
     * the cell before it had just written. The last cell has none after it, so
     * its insertion is laid once the row is done. */
    double            carried  = 0.0;
    scaled_row        above;
    weighing          weight;
    landing           at;
    double            confidence;
    hts_pos_t         shift;
    hts_pos_t         above_width;

    if (i == 0)
        return;

    above       = scaled_row_of(ctx, i - 1);
    weight      = weighing_of(ctx, i);
    at          = landing_of(ctx, i);
    confidence  = confidence_at(ctx, i);
    shift       = shift_between(ctx, i - 1, i);
    above_width = width_at(ctx, i - 1);

    for (hts_pos_t k = 0; k < width; k++) {
        hts_pos_t opening = k + shift;   /* this position, one row up */
        double    matched = forward_at(&front, k, STATE_MATCH);
        double    skipped = forward_at(&front, k, STATE_DELETION);
        double    paired  = matched * back[k][STATE_MATCH];
        double    passed  = skipped * back[k][STATE_DELETION];
        double    opened  = within(opening, above_width)
                          ? forward_at(&above, opening, STATE_MATCH)
                          * back[k][STATE_INSERTION]
                          : 0.0;

        at.coverage[k] += paired * confidence;
        at.spanned[k]  += paired + passed;

        /* What the cell before this one opened, and then what this one lays. */
        at.mutations[k] += carried;
        at.mutations[k] += weight.substitution * paired * terms[k].modification
                         + weight.deletion * skipped * pairings[k];

        carried = weight.insertion * opened;
    }

    at.mutations[width] += carried;
}

/* Forward and backward have to describe the same set of paths, or what they
 * multiply out to are not posteriors at all. Every row coming to one is exactly
 * that statement, and it is cheap enough to make of one of them. */
static bool normalized(const context *ctx)
{
    scaled_row       front = scaled_row_of(ctx, 0);
    const band_cell *back  = backward_row_of(ctx, 0);
    double           total = 0.0;

    for (hts_pos_t k = 0; k < width_at(ctx, 0); k++)
        for (int s = 0; s < N_STATES; s++)
            total += forward_at(&front, k, s) * back[k][s];

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
    free(scratch->terms);
    free(scratch->scale);
    free(scratch->backward);
    free(scratch->pairings);
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
    band_cell  *forward;
    band_cell  *backward;
    double     *pairings;
    cell_terms *terms;

    if (rows <= scratch->matrix_rows && widest <= scratch->widest)
        return 0;

    rows   = rows   > scratch->matrix_rows ? rows   : scratch->matrix_rows;
    widest = widest > scratch->widest      ? widest : scratch->widest;

    forward  = realloc(scratch->forward, rows * widest * sizeof *forward);
    backward = realloc(scratch->backward, 2 * widest * sizeof *backward);
    pairings = realloc(scratch->pairings, widest * sizeof *pairings);
    terms    = realloc(scratch->terms, rows * widest * sizeof *terms);

    if (forward)
        scratch->forward = forward;
    if (backward)
        scratch->backward = backward;
    if (pairings)
        scratch->pairings = pairings;
    if (terms)
        scratch->terms = terms;

    if (!forward || !backward || !pairings || !terms)
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

    for (size_t i = 0; i < ctx->rows; i++) {
        hts_pos_t width = width_at(ctx, i);

        if (width > widest)
            widest = width;
    }

    return widest;
}

/* Every position a cell of any row can name, and the one before the earliest of
 * them, that being where the first cell of a row lays what it pairs. A row's
 * reach is what the CIGAR path crosses on it give or take its own half-width,
 * and a row further along may be narrower than one behind it, so neither end
 * belongs to a particular row and both are looked for. */
static extent window_of(const context *ctx)
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

    return (extent){
        .origin = first - 1,
        .len    = (size_t)(last - first) + 1,
    };
}

/* The window is added to and never assigned, so what one read leaves in it must
 * not reach the next. */
static void clear_window(const context *ctx)
{
    phmm_scratch *scratch = ctx->scratch;
    size_t        len     = ctx->window.len;

    memset(scratch->coverage, 0, len * sizeof *scratch->coverage);
    memset(scratch->spanned, 0, len * sizeof *scratch->spanned);
    memset(scratch->mutations, 0, len * sizeof *scratch->mutations);
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

    ctx->window = window_of(ctx);

    if (grow_window(ctx->scratch, ctx->window.len) < 0)
        return false;

    clear_window(ctx);

    return true;
}

phmm_status phmm_run(const phmm *model, const phred *quality,
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

    if (!prepare(&ctx))
        return PHMM_NO_MEMORY;

    if (!forward(&ctx) || !backward(&ctx))
        return PHMM_UNSOUND;

    out->origin    = ctx.window.origin;
    out->len       = ctx.window.len;
    out->coverage  = scratch->coverage;
    out->spanned   = scratch->spanned;
    out->mutations = scratch->mutations;

    return PHMM_OK;
}
