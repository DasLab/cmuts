/* phmm.c -- a banded pair HMM over one read and the reference it was placed on.
 *
 * Three states: M pairs a read base with a reference base, I takes a read base
 * the reference has nothing for, D passes over a reference base the read has
 * nothing for. Both indices are prefix lengths, so cell (i, j) is the first i
 * placed read bases against the first j reference bases, and the state is the
 * last operation. A cell in M or I reads reference base j - 1, read base i - 1.
 *
 * Each read base carries only a band of reference positions: what the CIGAR
 * path crosses on that row, plus a half-width either side. Every probability
 * here is therefore conditioned on the alignment staying inside the band. That
 * is the approximation: local ambiguity is marginalized away, a grossly
 * misplaced read is not rescued. A band of 0 marginalizes over the CIGAR path
 * alone.
 *
 * The forward pass divides each row by its own total and the backward pass
 * divides by the same numbers, which lets a posterior be read off as the
 * product of the two with no separate normalizer and no logarithms. Those
 * divisors multiply back up to the read's likelihood, which is never formed.
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

/* Results of one cell's base comparison, used by four separate readers: the
 * forward pass, the backward pass, the deletion closing on the row, and the
 * accumulation. One comparison serves all four, so it is made once and stored.
 * Stored for the whole matrix, not one row, since the backward pass reads rows
 * the forward pass has already left behind. */
typedef struct {
    double emission;
    double modification;
} cell_terms;

/* The emission of a comparison that settles nothing: an inserted base, which
 * answers to no reference position, or a base neither side has named. */
#define UNINFORMATIVE (1.0 / NUC_BASES)

/* A misread or modified base lands on one of the bases the reference does not
 * name, with equal probability for each. Both this and UNINFORMATIVE derive
 * from the alphabet width, so nuc.h is the only place either changes. */
#define OTHER_BASES ((double)(NUC_BASES - 1))

/* Forward times backward must come to one on every row. A departure past this
 * indicates an index error, not rounding, and the read is discarded. */
#define NORMALIZATION_TOLERANCE 1e-6

/* Rows are stored at the widest row's stride, so a row is located by
 * multiplication and no offset table is needed. Only the loops are ragged, and
 * they are what the cost follows. */
struct phmm_scratch {
    aln_place *places;     /* one per placed read base, and one before them */
    double    *scale;      /* what each forward row is scaled by, which is the
                              reciprocal of what it came to */
    band_cell *forward;    /* rows, each of widest cells */
    cell_terms *terms;     /* one per cell of it */
    /* Only the current row and the one below it are read, so two suffice
     * however many rows the read has. */
    band_cell *backward;
    /* The value of a pairing on the row below, one per cell of the current
     * row. The backward pass forms it and the deletion closing there is a
     * product with it, so it is passed along instead of recomputed. */
    double    *pairings;
    double    *coverage;   /* the window handed back */
    double    *spanned;
    double    *mutations;
    size_t     rows;         /* rows places and scale are sized for */
    size_t     matrix_rows;  /* rows the forward matrix is sized for */
    size_t     widest;       /* cells a row of it holds */
    size_t     window;       /* positions the last three are sized for */
};

/* A stretch of the reference: where it begins and how far it runs. */
typedef struct {
    hts_pos_t origin;
    size_t    len;
} extent;

/* Everything one read is marginalized against.
 *
 * The caller's inputs come first and are fixed for the read. The rest is
 * derived from them and is undefined until prepare() has run; before that it
 * holds the previous read's band. prepare() is the only writer, which is why
 * every other function here takes a const context. */
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

/* A read base can disagree with the reference two ways -- the template was
 * modified, or the base was misread -- and only the second is what a quality
 * score measures. Agreement plus three disagreements sums to one, which makes
 * these a distribution over the four bases. */
static double agreement_chance(double modification, double error)
{
    return (1.0 - modification) * (1.0 - error)
         + modification * error / OTHER_BASES;
}

static double disagreement_chance(double modification, double error)
{
    return (1.0 - modification) * error / OTHER_BASES
         + modification * (1.0 - error / OTHER_BASES) / OTHER_BASES;
}

/* Given the base as read, the chance the template really differed here.
 *
 * A base agreeing with the reference may have been modified and misread back
 * into agreement; a disagreeing one may be an unmodified base misread. The
 * first is negligible, the second is not, which is why a poorly read
 * disagreement counts for less than a clean one. */
static double phmm_modification(const phmm *model, bool agree, double error)
{
    double modification = model->params.modification;

    return agree
         ? modification * error / OTHER_BASES
             / agreement_chance(modification, error)
         : modification * (1.0 - error / OTHER_BASES) / OTHER_BASES
             / disagreement_chance(modification, error);
}

/* ------------------------------------------------------------------------ */
/* One read base against the reference                                       */
/* ------------------------------------------------------------------------ */

/* A record storing no qualities weighs every base fully, leaving the
 * modification rate alone to explain a disagreement. */
static double error_at(const context *ctx, int32_t query)
{
    return ctx->read->qual
         ? phred_error(ctx->quality, ctx->read->qual[query])
         : 0.0;
}

/* The emission of an agreeing or disagreeing pairing, and how much of it is
 * attributable to a real difference in the template. */
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
 * A row pairs one read base, read and scored once, so a cell only selects among
 * three cases: the reference offers the same base, a different one, or nothing
 * comparable. All three are computed on entering the row, keeping the divisions
 * a modification costs out of the per-cell path.
 *
 * Where the read names no base, all three collapse to the same value. */
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

/* Cell (i, j) compares the row's base against reference base j - 1. A position
 * past either end of the reference is unknown, not a disagreement. */
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

/* Confidence in the base a row pairs, whatever it was set against: a base
 * neither side could name was still read, at the quality its score gives. */
static double confidence_at(const context *ctx, size_t i)
{
    return 1.0 - error_at(ctx, ctx->span.begin + (int32_t)i - 1);
}

/* ------------------------------------------------------------------------ */
/* The band                                                                  */
/* ------------------------------------------------------------------------ */

/* A forward row as stored, which is not the value it represents. See
 * scaled_row below before reading one. */
static band_cell *row_of(const context *ctx, size_t i)
{
    return ctx->scratch->forward + i * (size_t)ctx->widest;
}

/* The same row, for readers.
 *
 * A pointer to an array takes its own qualifier only from C23, so a read-only
 * row cannot simply be declared const at each use. The cast lives here and in
 * read_backward_row_of instead of at every call site. */
static const band_cell *read_row_of(const context *ctx, size_t i)
{
    return (const band_cell *)row_of(ctx, i);
}

/* A forward row together with the factor that restores it.
 *
 * A row is left scaled by every divisor except its own, that last one being
 * folded into the transitions leaving the row instead of walked over the row a
 * second time. What is stored is therefore the row's total times its true
 * value; the two are a probability only together.
 *
 * Undo it here and nowhere else. Reading a forward row directly is the one
 * error this file cannot detect: the numbers look like posteriors but are wrong
 * by a per-row factor, so nothing overflows, nothing fails to sum, and the
 * counts are silently wrong.
 *
 * The row and its divisor are fetched together, once. A pass that wrote
 * anything between the two fetches would have to repeat both: the window it
 * writes is doubles, as is a row, and the compiler cannot assume they do not
 * alias. */
typedef struct {
    const band_cell *cell;
    double           scale;
} scaled_row;

static scaled_row scaled_row_of(const context *ctx, size_t i)
{
    return (scaled_row){
        .cell  = read_row_of(ctx, i),
        .scale = ctx->scratch->scale[i],
    };
}

/* The true value of one cell. */
static double forward_at(const scaled_row *row, hts_pos_t k, int state)
{
    return row->cell[k][state] * row->scale;
}

static cell_terms *terms_of(const context *ctx, size_t i)
{
    return ctx->scratch->terms + i * (size_t)ctx->widest;
}

/* Row parity selects which of the two buffers it occupies. */
static band_cell *backward_row_of(const context *ctx, size_t i)
{
    return ctx->scratch->backward + (i & 1) * (size_t)ctx->widest;
}

/* As read_row_of, for the backward buffers. */
static const band_cell *read_backward_row_of(const context *ctx, size_t i)
{
    return (const band_cell *)backward_row_of(ctx, i);
}

/* Reference the row's base accounts for beyond the one it pairs with, which is
 * whatever a deletion following it passes over. */
static hts_pos_t skip_at(const context *ctx, size_t i)
{
    return ctx->scratch->places[i].last - ctx->scratch->places[i].first;
}

/* A row covers everything the CIGAR path crosses on it, plus a half-width
 * either side -- not a half-width either side of a single point. Deletions are
 * why: the path runs along a row to cross one, so a row with a deletion after
 * it is a stretch that even a band of 0 holds whole. A read is never denied the
 * alignment it arrived with; the band only adds room to depart from it. */
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
 * follow the CIGAR and may widen or narrow, so the offset between two of them
 * is the distance between their first cells. */
static hts_pos_t shift_between(const context *ctx, size_t from, size_t to)
{
    return origin_of(ctx, to) - origin_of(ctx, from);
}

/* Rows differ in width, so a cell index is meaningful only against its own
 * row's width. Checking it against another row's is how the two passes would
 * come to disagree about which paths exist.
 *
 * The width is a parameter because a row's width is fixed while it is walked,
 * and re-fetching it per cell is most of what this test would cost. */
static bool within(hts_pos_t k, hts_pos_t width)
{
    return k >= 0 && k < width;
}

/* ------------------------------------------------------------------------ */
/* Forward                                                                   */
/* ------------------------------------------------------------------------ */

/* Nothing has been read yet, so each cell of the first row is the alignment
 * about to begin at that reference position, all equally likely. None begins
 * with a deletion: reference passed over before the read starts carries no
 * information and its position could never be recovered. */
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
 * everything constant along the row multiplied out once: the row above's scale
 * factor, and the flat emission an inserted base carries. Only the pairing's
 * own emission is left for the cell to supply.
 *
 * Folding the scaling in here saves a second walk over the row. Dividing the
 * row by its total and stepping from it undivided give the same numbers, but
 * that costs a multiply, load and store per state per cell, against a handful
 * of multiplies for the whole row here.
 *
 * This is the only place that reads a forward row without forward_at, and it
 * performs the same correction by another route: the divisor that would
 * multiply every cell above is carried instead by these five constants. See
 * forward_at for why the correction is needed. */
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

/* The two states that consume a read base, and so read from the row above. */
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

/* A deletion advances along the reference without the read moving, so it is the
 * only state whose row depends on itself: a cell reads the cell to its left,
 * which must already hold this row's value, not the previous row's. */
static double deleted_from(const phmm *model, double left_match,
                           double left_deletion)
{
    return model->match_to_deletion    * left_match
         + model->deletion_to_deletion * left_deletion;
}

/* A deletion may neither open nor close a read: reference passed over before
 * the first read base or after the last has no recoverable position, and
 * allowing it would attribute an event to nothing in the read.
 *
 * Both passes call this, rather than each enforcing the rule separately. A step
 * one allows and the other forbids is a disagreement about which paths exist,
 * and posteriors read off the two would not sum to one. */
static bool deletions_live(const context *ctx, size_t i)
{
    return i > 0 && i + 1 < ctx->rows;
}

/* One row of the forward pass, returning its unscaled total.
 *
 * Every input a cell needs is available when it is reached: its comparison is
 * made here, the two read-consuming states read the row above, and the deletion
 * reads the cell to the left, written one iteration ago. That makes this a
 * single pass over the row instead of four, which at a narrow band is most of
 * what a row costs. The total accumulates in cell order, so it is the same sum
 * the scaling afterwards divides by. */
static double forward_row(const context *ctx, size_t i)
{
    const phmm      *model = ctx->model;
    descent          step  = descent_into(ctx, i);
    band_cell       *row   = row_of(ctx, i);
    const band_cell *above = read_row_of(ctx, i - 1);
    cell_terms      *terms = terms_of(ctx, i);
    row_terms        each  = row_terms_of(ctx, i);
    hts_pos_t        shift = shift_between(ctx, i - 1, i);
    hts_pos_t        width = width_at(ctx, i);
    hts_pos_t        above_width = width_at(ctx, i - 1);
    bool             live  = deletions_live(ctx, i);
    /* The cell to the left, held in registers instead of re-read from the row
     * just written: the deletion chain runs the width of the row, and a load
     * waiting on the preceding store lengthens every link of it. */
    double           left_match    = 0.0;
    double           left_deletion = 0.0;
    /* One running sum per state, making three addition chains as deep as the
     * row is wide instead of one chain three times that.
     *
     * This is the only place here that sums a row out of cell order. The result
     * is the row's divisor, and nothing downstream depends on the order, though
     * it is not bit-exact in principle. */
    double           total_paired   = 0.0;
    double           total_inserted = 0.0;
    double           total_deleted  = 0.0;

    for (hts_pos_t k = 0; k < width; k++) {
        hts_pos_t diagonal = k - 1 + shift;  /* a position back, one row up */
        hts_pos_t straight = k + shift;      /* this position, one row up */
        double    paired, inserted, deleted;

        terms[k] = terms_at(ctx, &each, position_of(ctx, i, k));

        paired = within(diagonal, above_width)
               ? paired_from(&step, above[diagonal], terms[k].emission)
               : 0.0;
        inserted = within(straight, above_width)
                 ? inserted_from(&step, above[straight])
                 : 0.0;
        deleted = live && k > 0
                ? deleted_from(model, left_match, left_deletion)
                : 0.0;

        row[k][STATE_MATCH]     = paired;
        row[k][STATE_INSERTION] = inserted;
        row[k][STATE_DELETION]  = deleted;

        left_match    = paired;
        left_deletion = deleted;

        total_paired   += paired;
        total_inserted += inserted;
        total_deleted  += deleted;
    }

    return (total_paired + total_inserted) + total_deleted;
}

/* Every row is scaled by its own total, keeping the numbers near one however
 * long the read: unscaled, a forward pass underflows a double within a few
 * hundred bases. The scaling itself happens in the descent out of the row, so
 * only the factor it needs is stored here. */
static bool record_total(const context *ctx, size_t i, double total)
{
    if (!(total > 0.0) || !isfinite(total))
        return false;

    ctx->scratch->scale[i] = 1.0 / total;

    return true;
}

/* The first row is initialized directly, having no row above to step from. */
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
 * without passing over a reference base. There is no row below for a deletion
 * to close on, so that state is zero throughout. */
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

/* Everything reachable from a cell. For the two states that stay on this row,
 * that is the cell to its right. Only transitions crossing to the row below
 * carry that row's scaling; the deletion chain along this row was scaled as it
 * was written.
 *
 * The value of a pairing on the row below is stored for the accumulation, the
 * deletion closing there being a product with exactly that. */
static void backward_row(const context *ctx, size_t i)
{
    const phmm *model = ctx->model;
    hts_pos_t   shift = shift_between(ctx, i, i + 1);
    double      below_scale = ctx->scratch->scale[i + 1];
    bool        live  = deletions_live(ctx, i);

    band_cell        *row      = backward_row_of(ctx, i);
    const band_cell  *below    = read_backward_row_of(ctx, i + 1);
    const cell_terms *terms    = terms_of(ctx, i + 1);
    double           *pairings = ctx->scratch->pairings;
    hts_pos_t         width    = width_at(ctx, i);
    hts_pos_t         below_width = width_at(ctx, i + 1);
    /* The cell to the right, held for the same reason the forward pass holds
     * the cell to its left. */
    double            right_deletion = 0.0;

    for (hts_pos_t k = width; k-- > 0; ) {
        hts_pos_t diagonal = k + 1 - shift;  /* a position on, one row down */
        hts_pos_t straight = k - shift;      /* this position, one row down */
        double    paired   = 0.0;
        double    inserted = 0.0;
        double    deleted  = live && k + 1 < width ? right_deletion : 0.0;

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

        right_deletion = row[k][STATE_DELETION];
    }
}

/* ------------------------------------------------------------------------ */
/* What a row says about the reference                                       */
/* ------------------------------------------------------------------------ */

/* The three window fields, each advanced to where a row's first cell enters
 * the window, so a cell is addressed by its own index with no bounds check.
 *
 * A cell at k denotes reference prefix length j: it pairs the base before it
 * and lays an insertion at j, so a row touches the window over
 * [j0 - 1, j0 + width - 1] and nowhere else. window_of sizes the window from
 * one before the earliest position any row names to the furthest, so every
 * such position is inside it. That a row is a contiguous stretch of the window
 * at a fixed offset is an invariant of window_of, not of this function.
 *
 * Positions past either end of the reference are kept. Some of what the band
 * names really does lie outside, since a read near a boundary has paths that
 * leave it; the caller decides which to keep, and deciding it here as well
 * would leave two places to change. */
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

/* What one event of each kind contributes to the position it is laid at, per
 * unit of posterior.
 *
 * Each is the kind's weight times every factor of the event's posterior that
 * does not vary by cell: for a deletion, the step back out of the run; for an
 * insertion, the step into it and the flat emission it carries.
 *
 * An insertion also carries the divisor of the row it opens on, as every
 * posterior over a step between two rows does. That comes from the recursion,
 * not from how a row is stored -- the latter is undone by forward_at, and the
 * two are easily confused.
 *
 * Valid only for a row with one above it; an insertion on the first row is one
 * nothing carried into. */
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

/* How each event contributes to the three per-position quantities.
 *
 * A pairing spans the position it pairs with, covers it in proportion to the
 * confidence in the base, and contributes the part of its posterior belonging
 * to a real difference in the template. A deletion spans every position it
 * passes over but covers none, there being no base read there; it counts as a
 * modification only where it opened, one adduct stopping one reverse
 * transcriptase once whatever length it then skipped. An insertion counts the
 * same way, at the position it precedes, and neither spans nor covers.
 *
 * Only the mutations are weighted by event kind. How a position was reached
 * does not bear on whether it was reached.
 *
 * A deletion is laid at the end of its run, not its start, because reverse
 * transcription reads the template from the 3' end: the last base a deletion
 * passes over is the first the enzyme met and the likeliest to carry what
 * stopped it. Both are marginals of one joint over runs, so the choice moves
 * where a deletion is counted, not how much of it there is. A deletion ends
 * against the pairing on the row below, left by the backward pass; an insertion
 * opens out of the pairing on the row above, still where the forward pass
 * wrote it.
 *
 * The row's scaling is undone on the two states read here, which is cheaper
 * than undoing it across the whole row as it was written.
 *
 * The first row is the alignment about to begin. Its posterior says where the
 * read starts, not what any base was set against, so it contributes nothing. */
static void accumulate_row(const context *ctx, size_t i)
{
    phmm_scratch     *scratch  = ctx->scratch;
    scaled_row        front    = scaled_row_of(ctx, i);
    const band_cell  *back     = read_backward_row_of(ctx, i);
    const cell_terms *terms    = terms_of(ctx, i);
    const double     *pairings = scratch->pairings;
    hts_pos_t         width    = width_at(ctx, i);
    /* An insertion belongs at the position its cell precedes, which is where
     * the next cell lays what it pairs. Carrying it to the following iteration
     * keeps each cell writing one mutations slot and no other; laying it ahead
     * would have every cell read back a slot the previous one just wrote. The
     * last cell has no successor, so its insertion is laid after the loop. */
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

        /* What the previous cell opened, then what this one lays. */
        at.mutations[k] += carried;
        at.mutations[k] += weight.substitution * paired * terms[k].modification
                         + weight.deletion * skipped * pairings[k];

        carried = weight.insertion * opened;
    }

    at.mutations[width] += carried;
}

/* Forward and backward must describe the same set of paths, or their product
 * is not a posterior. Every row summing to one is that statement, and checking
 * a single row is cheap. */
static bool normalized(const context *ctx)
{
    scaled_row       front = scaled_row_of(ctx, 0);
    const band_cell *back  = read_backward_row_of(ctx, 0);
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

/* A buffer that grew is kept even if a later one fails, so a scratch that fails
 * to grow is still usable, just larger in places than it needs to be. Nothing
 * shrinks; a scratch is reused across every read one worker sees.
 *
 * The row buffers grow first and alone: their count comes from the CIGAR, and
 * the places must be written before the widest row is known. */
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

/* Sized by rows times widest row, so a band wider than the previous read's
 * triggers a regrow however few rows this one needs. */
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
 * the buffers must accommodate. At least one, so an absurd band still leaves an
 * addressable matrix; its rows are then empty and the pass fails on the
 * first. */
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

/* Every position any cell can name, plus the one before the earliest, where a
 * row's first cell lays what it pairs. A row reaches what the CIGAR path
 * crosses on it, give or take its half-width, and a later row may be narrower
 * than an earlier one, so neither end belongs to a fixed row and both are
 * searched for. */
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

/* The window is accumulated into, never assigned, so it must be cleared
 * between reads. */
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
