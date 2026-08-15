/* phmm.c -- a banded pair HMM over one read and the reference it was placed on.
 *
 * Three states: M pairs a read base with a reference base, I consumes a read
 * base with no reference base, D consumes a reference base with no read base.
 * Both indices are prefix lengths, so cell (i, j) is the first i placed read
 * bases against the first j reference bases, with the state giving the last
 * operation. A cell in M or I compares reference base j - 1 with read base
 * i - 1.
 *
 * Each read base carries only a band of reference positions: what the CIGAR
 * path crosses on that row, plus a half-width either side. Every probability
 * here is conditioned on the alignment staying inside the band, so local
 * ambiguity is marginalized away but a grossly misplaced read is not corrected.
 * A band of 0 marginalizes over the CIGAR path only.
 *
 * Both passes divide each row by the forward row's total, so a posterior is the
 * product of the two with no separate normalizer and no logarithms.
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

/* One cell of a band row. A named type so that a row of cells can be returned. */
typedef double band_cell[N_STATES];

/* One cell's base comparison, read by the forward pass, the backward pass and
 * the accumulation. Stored for the whole matrix, since the backward pass reads
 * rows the forward pass has already gone past. */
typedef struct {
    double emission;
    double modification;
} cell_terms;

/* Emission of a comparison that distinguishes nothing: an inserted base, or a
 * base one side does not hold. */
#define UNINFORMATIVE (1.0 / NUC_BASES)

/* A misread or modified base is equally likely to be any base the reference does
 * not hold. Derived from the alphabet width, as UNINFORMATIVE is. */
#define OTHER_BASES ((double)(NUC_BASES - 1))

/* Forward times backward sums to one on the first row. A finite departure past this
 * indicates an index error and not rounding, and the run stops. */
#define NORMALIZATION_TOLERANCE 1e-6

/* Rows are stored at the widest row's stride, so a row is located by
 * multiplication. Only the loops are ragged. */
struct phmm_scratch {
    aln_place *places;     /* one per placed read base, and one before them */
    double    *scale;      /* reciprocal of each forward row's total */
    band_cell *forward;    /* rows, each of widest cells */
    cell_terms *terms;     /* one per cell of it */
    /* Two rows suffice: only the current row and the one below it are read. */
    band_cell *backward;
    double    *coverage;   /* the window returned to the caller */
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
 * derived from them by prepare(), which is the only writer, so every other
 * function here takes a const context. */
typedef struct {
    const phmm            *model;
    const phred           *quality;
    const cm_bam_record   *read;
    const cm_fasta_record *ref;
    phmm_scratch          *scratch;
    const int             *half;    /* half-width either side of the CIGAR, one
                                       per row; the caller's */

    aln_span               span;    /* the stretch of the read that is placed */
    size_t                 rows;    /* placed bases, and one row before them */
    hts_pos_t              widest;  /* cells the widest row holds */
    extent                 window;  /* every position any row writes to */
} context;

/* ------------------------------------------------------------------------ */
/* The model                                                                 */
/* ------------------------------------------------------------------------ */

phmm_weights phmm_default_weights(void)
{
    return (phmm_weights){
        .weight = {
            [PHMM_SUBSTITUTION] = 1.0,
            [PHMM_DELETION]     = 1.0,
            [PHMM_INSERTION]    = 0.0,
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

/* Returns the chance a read base agrees with the reference base it was templated from,
 * and the chance it differs. A base can disagree two ways -- the template was
 * modified, or the base was misread -- and a quality score measures only the
 * second. Agreement plus three disagreements sums to one. */
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

/* Gives the chance the template differed here, given the base as read. An agreeing base
 * may have been modified and misread back into agreement; a disagreeing one may
 * be an unmodified base misread, so a poorly read disagreement counts for less
 * than a clean one. */
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

/* Returns the chance the base at the given query offset was misread. A record storing
 * no qualities weighs every base fully, leaving the modification rate to explain a
 * disagreement. */
static double error_at(const context *ctx, int32_t query)
{
    return ctx->read->qual
         ? phred_error(ctx->quality, ctx->read->qual[query])
         : 0.0;
}

/* Gives the emission of an agreeing or disagreeing pairing, and how much of it is
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

/* The three values a cell at a reference position can take: the reference
 * holds the same base, a different one, or nothing comparable. All three are
 * computed on entering the row, keeping the divisions in phmm_modification out
 * of the per-cell path. Where the read holds no base they are equal. */
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

/* Returns the comparison for a cell at reference prefix length j, which pairs the row's
 * base with reference base j - 1. A position past either end of the reference has
 * no base to pair with, so its emission is zero; an ambiguous base is a real base
 * of unknown identity and takes the uninformative emission. */
static cell_terms terms_at(const context *ctx, const row_terms *row,
                           hts_pos_t j)
{
    nuc theirs;

    if (j < 1 || (size_t)j > ctx->ref->len) {
        return (cell_terms){ .emission = 0.0, .modification = 0.0 };
    }

    theirs = nuc_from_char(ctx->ref->seq[j - 1]);

    if (!nuc_is_base(theirs)) {
        return row->neither;
    }

    return theirs == row->ours ? row->agree : row->differ;
}

/* Confidence in the base a row pairs, regardless of what it was compared against: an
 * unnamed base was still read, at the quality its score gives. */
static double confidence_at(const context *ctx, size_t i)
{
    return 1.0 - error_at(ctx, ctx->span.begin + (int32_t)i - 1);
}

/* ------------------------------------------------------------------------ */
/* The band                                                                  */
/* ------------------------------------------------------------------------ */

/* Returns a forward row as stored, which is not the value it represents. See
 * scaled_row. */
static band_cell *row_of(const context *ctx, size_t i)
{
    return ctx->scratch->forward + i * (size_t)ctx->widest;
}

/* Returns the same row, const. A pointer to an array takes its own qualifier only from
 * C23, so the cast lives here and in read_backward_row_of, and not at every
 * call site. */
static const band_cell *read_row_of(const context *ctx, size_t i)
{
    return (const band_cell *)row_of(ctx, i);
}

/* A forward row together with the factor that restores it.
 *
 * A row is stored scaled by every divisor except its own, that last one being
 * folded into the transitions leaving the row. What is stored is therefore the
 * row's total times its true value, so the two are a probability only together.
 *
 * Reading a forward row without undoing the scaling is the one error this file
 * cannot detect: the values look like posteriors but are wrong by a per-row
 * factor, so nothing overflows and no row fails to sum. */
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

/* Returns the true value of one cell. */
static double forward_at(const scaled_row *row, hts_pos_t k, int state)
{
    return row->cell[k][state] * row->scale;
}

static cell_terms *terms_of(const context *ctx, size_t i)
{
    return ctx->scratch->terms + i * (size_t)ctx->widest;
}

/* Returns the buffer holding row i, chosen by its parity. */
static band_cell *backward_row_of(const context *ctx, size_t i)
{
    return ctx->scratch->backward + (i & 1) * (size_t)ctx->widest;
}

/* As read_row_of, for the backward buffers. */
static const band_cell *read_backward_row_of(const context *ctx, size_t i)
{
    return (const band_cell *)backward_row_of(ctx, i);
}

/* Reference bases the row's base accounts for beyond the one it pairs with,
 * which is whatever a following deletion passes over. */
static hts_pos_t skip_at(const context *ctx, size_t i)
{
    return ctx->scratch->places[i].last - ctx->scratch->places[i].first;
}

/* Returns the cells in row i: the whole stretch the CIGAR path crosses on it, plus a
 * half-width either side. A row with a deletion after it spans the whole skip
 * even at a band of 0, so the band only adds room to depart from the reported
 * alignment. */
static hts_pos_t width_at(const context *ctx, size_t i)
{
    return skip_at(ctx, i) + 2 * (hts_pos_t)ctx->half[i] + 1;
}

/* Returns the reference prefix length a row's first cell stands for. */
static hts_pos_t origin_of(const context *ctx, size_t i)
{
    return ctx->scratch->places[i].first - ctx->half[i];
}

/* Returns the reference prefix length a cell of a row stands for. */
static hts_pos_t position_of(const context *ctx, size_t i, hts_pos_t k)
{
    return origin_of(ctx, i) + k;
}

/* Returns what to add to a cell index to reach the same position on another row. Rows
 * may widen or narrow, so the offset is the distance between their first cells. */
static hts_pos_t shift_between(const context *ctx, size_t from, size_t to)
{
    return origin_of(ctx, to) - origin_of(ctx, from);
}

/* Returns whether k is a cell of a row of the given width. Rows differ in width, so an
 * index is meaningful only against its own row's. The width is a parameter,
 * since it is fixed while a row is walked. */
static bool within(hts_pos_t k, hts_pos_t width)
{
    return k >= 0 && k < width;
}

/* ------------------------------------------------------------------------ */
/* Forward                                                                   */
/* ------------------------------------------------------------------------ */

/* Fills the first row and returns its total. Each cell is the alignment beginning
 * at that reference position, all equally likely; none begins with a deletion,
 * reference passed over before the read starts having no recoverable position. */
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
 * factor, and the flat emission an inserted base carries. The cell supplies only
 * the pairing's own emission. Folding the scaling in here saves a second walk
 * over the row.
 *
 * This is the only place a forward row is read without forward_at. The same
 * correction is applied by another route: the divisor that would multiply every
 * cell above is carried by these five constants instead. */
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

/* Returns a cell's two states that consume a read base, and so step from the row above. */
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

/* Returns a cell's deletion state, stepped from the cell to its left. A deletion
 * advances along the reference without the read moving, so this is the only state that
 * depends on its own row. */
static double deleted_from(const phmm *model, double left_match,
                           double left_deletion)
{
    return model->match_to_deletion    * left_match
         + model->deletion_to_deletion * left_deletion;
}

/* Returns whether row i may hold a deletion. The first and last rows may not: reference
 * passed over before the first read base or after the last has no recoverable
 * position. */
static bool deletions_live(const context *ctx, size_t i)
{
    return i > 0 && i + 1 < ctx->rows;
}

/* Fills row i and returns its unscaled total. Every input a cell needs is
 * available when it is reached -- its comparison is made here, the two
 * read-consuming states read the row above, and the deletion reads the cell to
 * the left, written one iteration earlier -- so one pass over the row suffices
 * where four would otherwise be needed. */
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
    /* The cell to the left, held in a local: the deletion chain runs the width of
     * the row, and a load waiting on the preceding store would lengthen every link
     * of it. */
    double           left_match    = 0.0;
    double           left_deletion = 0.0;
    /* One running sum per state, giving three addition chains the width of the
     * row in place of one chain three times as deep. */
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

/* Stores row i's scale factor. Every row is scaled by its own total, an unscaled forward
 * pass underflowing a double within a few hundred bases; the scaling itself is applied in
 * the descent out of the row.
 *
 * A total that is zero or subnormal has lost the precision the scaling needs: the band
 * admits no path this read can be scored along, which is the read against the rates and
 * not a fault. A total that is negative or not a number is one. */
static phmm_status record_total(const context *ctx, size_t i, double total)
{
    if (total < 0.0 || !isfinite(total)) {
        return PHMM_UNSOUND;
    }

    if (!isnormal(total)) {
        return PHMM_NO_PATH;
    }

    ctx->scratch->scale[i] = 1.0 / total;

    return PHMM_OK;
}

/* Runs the forward pass. The first row is filled directly, having no row above to step
 * from. */
static phmm_status forward(const context *ctx)
{
    phmm_status status = record_total(ctx, 0, forward_first_row(ctx));

    for (size_t i = 1; status == PHMM_OK && i < ctx->rows; i++) {
        status = record_total(ctx, i, forward_row(ctx, i));
    }

    return status;
}

/* ------------------------------------------------------------------------ */
/* Accumulating a row into the window                                        */
/* ------------------------------------------------------------------------ */

/* The three window fields, each advanced to where a row's first cell enters the
 * window, so a cell is addressed by its own index with no bounds check.
 *
 * A cell at k denotes reference prefix length j: it pairs the base before it and
 * lays an insertion at j, so a row touches the window over
 * [j0 - 1, j0 + width - 1]. That every such position lies inside the window,
 * contiguously and at a fixed offset, is an invariant of window_of.
 *
 * The band is not clamped, so a row near either end of the reference covers a few
 * positions outside it. Those are addressed like any other and never written
 * to. */
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

/* What one event of each kind contributes to the position it is laid at, per unit
 * of posterior: the kind's weight times every factor of the event's posterior
 * that does not vary by cell. For a deletion that is the step back out of the
 * run; for an insertion, the step into it and the flat emission it carries.
 *
 * An insertion also carries the divisor of the row it opens on, as every
 * posterior over a step between two rows does. That comes from the recursion and
 * not from how a row is stored, the latter being undone by forward_at.
 *
 * Valid only for a row with one above it, no insertion opening on the first
 * row. */
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

/* Everything the accumulation of one row works from, and the one position of it
 * still open.
 *
 * Cells arrive right to left, the order the backward pass forms them in. An
 * insertion belongs at the position the next cell to the right pairs, so that
 * position is complete only once the cell to its left is reached: what the cell
 * holds of its own is kept here for exactly one step, and written once whatever
 * carries into it is known. Each cell then writes one slot of each field and
 * reads none back.
 *
 * The held values start at zero, the rightmost cell having no cell to its right.
 *
 * Valid only for a row with one above it. */
typedef struct {
    scaled_row        front;
    scaled_row        above;
    const cell_terms *terms;
    weighing          weight;
    landing           at;
    double            confidence;
    hts_pos_t         up;           /* to the same position, one row up */
    hts_pos_t         above_width;
    /* Held from the cell to the right, without the insertion carried into its
     * position. */
    double            coverage;
    double            spanned;
    double            mutations;
} accumulation;

static accumulation accumulation_of(const context *ctx, size_t i)
{
    return (accumulation){
        .front       = scaled_row_of(ctx, i),
        .above       = scaled_row_of(ctx, i - 1),
        .terms       = terms_of(ctx, i),
        .weight      = weighing_of(ctx, i),
        .at          = landing_of(ctx, i),
        .confidence  = confidence_at(ctx, i),
        .up          = shift_between(ctx, i - 1, i),
        .above_width = width_at(ctx, i - 1),
    };
}

/* Writes the position completed by reaching cell k, and holds what this cell
 * contributes to the position to its left.
 *
 * A pairing spans the position it pairs with, covers it in proportion to the
 * confidence in the base, and contributes the part of its posterior belonging to
 * a real difference in the template. A deletion spans every position it passes
 * over and covers none, since no base was read there; it counts as a
 * modification only where it opened, one adduct having stopped one reverse
 * transcriptase whatever length was then skipped. An insertion counts the same
 * way, at the position it precedes, and covers nothing.
 *
 * A pairing and a deletion span their posterior unweighted. An insertion spans
 * what its weight gives, the same quantity it lays as a mutation: how far an
 * inserted base bears on a reference position is exactly how far it is taken to
 * be a modification there. Spanning it unweighted would make a weight of zero
 * into evidence against a modification, and not the absence of evidence.
 *
 * A deletion is laid at the end of its run and not its start, since reverse
 * transcription reads the template from the 3' end: the last base a deletion
 * passes over is the first the enzyme met. Both are marginals of one joint over
 * runs, so the choice moves where a deletion is counted, not how much there is.
 * A deletion closes against the pairing on the row below, which is why that is
 * passed in; an insertion opens out of the pairing on the row above.
 *
 * The row's scaling is undone on the two states read here, and not across the
 * whole row as it was written. */
static void accumulate_cell(accumulation *acc, hts_pos_t k,
                            const double *back, double pairing)
{
    hts_pos_t opening = k + acc->up;   /* this position, one row up */
    double    matched = forward_at(&acc->front, k, STATE_MATCH);
    double    skipped = forward_at(&acc->front, k, STATE_DELETION);
    double    paired  = matched * back[STATE_MATCH];
    double    passed  = skipped * back[STATE_DELETION];
    double    carried = within(opening, acc->above_width)
                      ? acc->weight.insertion
                      * forward_at(&acc->above, opening, STATE_MATCH)
                      * back[STATE_INSERTION]
                      : 0.0;

    acc->at.coverage[k + 1]  += acc->coverage;
    acc->at.spanned[k + 1]   += acc->spanned + carried;
    acc->at.mutations[k + 1] += acc->mutations + carried;

    acc->coverage  = paired * acc->confidence;
    acc->spanned   = paired + passed;
    acc->mutations = acc->weight.substitution * paired
                   * acc->terms[k].modification
                   + acc->weight.deletion * skipped * pairing;
}

/* Writes the last position held. The leftmost cell has nothing to its left, so
 * what it holds is the whole of that position. */
static void accumulate_end(const accumulation *acc)
{
    acc->at.coverage[0]  += acc->coverage;
    acc->at.spanned[0]   += acc->spanned;
    acc->at.mutations[0] += acc->mutations;
}

/* ------------------------------------------------------------------------ */
/* Backward                                                                  */
/* ------------------------------------------------------------------------ */

/* Fills the last row of the backward pass, accumulating into the window as it is formed
 * so that a cell's values stay in registers and are not read back.
 *
 * The alignment ends on this row, having paired or inserted its final base
 * without passing over a reference base. There is no row below for a deletion to
 * close on, so that state is zero throughout. */
static void backward_last_row(const context *ctx)
{
    size_t       i   = ctx->rows - 1;
    band_cell   *row = backward_row_of(ctx, i);
    accumulation acc = accumulation_of(ctx, i);
    double       cell[N_STATES] = {
        [STATE_MATCH]     = 1.0,
        [STATE_INSERTION] = 1.0,
        [STATE_DELETION]  = 0.0,
    };

    for (hts_pos_t k = width_at(ctx, i); k-- > 0; ) {
        row[k][STATE_MATCH]     = cell[STATE_MATCH];
        row[k][STATE_INSERTION] = cell[STATE_INSERTION];
        row[k][STATE_DELETION]  = cell[STATE_DELETION];

        accumulate_cell(&acc, k, cell, 0.0);
    }

    accumulate_end(&acc);
}

/* Fills row i of the backward pass and accumulates it. A cell holds everything
 * reachable from it, which for the two states staying on this row is the cell to
 * its right. Only transitions into the row below carry that row's scaling; the
 * deletion chain along this row was scaled as it was written.
 *
 * Valid only for a row with one above and one below, which is every row a deletion
 * may live on, so the deletion state is unconditional here. */
static void backward_row(const context *ctx, size_t i)
{
    const phmm *model = ctx->model;

    band_cell        *row         = backward_row_of(ctx, i);
    const band_cell  *below       = read_backward_row_of(ctx, i + 1);
    const cell_terms *below_terms = terms_of(ctx, i + 1);
    double            below_scale = ctx->scratch->scale[i + 1];
    hts_pos_t         below_width = width_at(ctx, i + 1);
    accumulation      acc         = accumulation_of(ctx, i);
    hts_pos_t         shift       = shift_between(ctx, i, i + 1);
    hts_pos_t         width       = width_at(ctx, i);
    /* The cell to the right, held for the same reason the forward pass holds
     * the cell to its left. */
    double            right_deletion = 0.0;

    for (hts_pos_t k = width; k-- > 0; ) {
        hts_pos_t diagonal = k + 1 - shift;  /* a position on, one row down */
        hts_pos_t straight = k - shift;      /* this position, one row down */
        double    pairing  = 0.0;
        double    inserted = 0.0;
        double    deleted  = k + 1 < width ? right_deletion : 0.0;
        double    cell[N_STATES];

        if (within(diagonal, below_width)) {
            pairing = below_terms[diagonal].emission
                    * below[diagonal][STATE_MATCH] * below_scale;
        }

        if (within(straight, below_width)) {
            inserted = UNINFORMATIVE
                     * below[straight][STATE_INSERTION] * below_scale;
        }

        cell[STATE_MATCH] = model->match_to_match     * pairing
                          + model->match_to_insertion * inserted
                          + model->match_to_deletion  * deleted;

        cell[STATE_INSERTION] = model->insertion_to_match     * pairing
                              + model->insertion_to_insertion * inserted;

        cell[STATE_DELETION] = model->deletion_to_match    * pairing
                             + model->deletion_to_deletion * deleted;

        row[k][STATE_MATCH]     = cell[STATE_MATCH];
        row[k][STATE_INSERTION] = cell[STATE_INSERTION];
        row[k][STATE_DELETION]  = cell[STATE_DELETION];

        accumulate_cell(&acc, k, cell, pairing);

        right_deletion = cell[STATE_DELETION];
    }

    accumulate_end(&acc);
}

/* Fills the first row of the backward pass, which is the alignment about to begin. Its
 * posterior gives where the read starts and not what any base was compared
 * against, so it contributes nothing and is written only for normalized() to
 * check. No deletion opens before the read starts, so that state is zero
 * throughout. */
static void backward_first_row(const context *ctx)
{
    const phmm *model = ctx->model;

    band_cell        *row         = backward_row_of(ctx, 0);
    const band_cell  *below       = read_backward_row_of(ctx, 1);
    const cell_terms *below_terms = terms_of(ctx, 1);
    double            below_scale = ctx->scratch->scale[1];
    hts_pos_t         below_width = width_at(ctx, 1);
    hts_pos_t         shift       = shift_between(ctx, 0, 1);
    hts_pos_t         width       = width_at(ctx, 0);

    for (hts_pos_t k = 0; k < width; k++) {
        hts_pos_t diagonal = k + 1 - shift;
        hts_pos_t straight = k - shift;
        double    pairing  = 0.0;
        double    inserted = 0.0;

        if (within(diagonal, below_width)) {
            pairing = below_terms[diagonal].emission
                    * below[diagonal][STATE_MATCH] * below_scale;
        }

        if (within(straight, below_width)) {
            inserted = UNINFORMATIVE
                     * below[straight][STATE_INSERTION] * below_scale;
        }

        row[k][STATE_MATCH] = model->match_to_match     * pairing
                            + model->match_to_insertion * inserted;

        row[k][STATE_INSERTION] = model->insertion_to_match     * pairing
                                + model->insertion_to_insertion * inserted;

        row[k][STATE_DELETION] = 0.0;
    }
}

/* Returns the first row's forward times backward, which is one where the two passes
 * describe the same set of paths.
 *
 * The first row is the last the backward pass computes, so it carries the whole backward
 * chain and every forward scale factor. It is also one of the two rows the pass keeps. */
static double agreement(const context *ctx)
{
    scaled_row       front = scaled_row_of(ctx, 0);
    const band_cell *back  = read_backward_row_of(ctx, 0);
    double           total = 0.0;

    for (hts_pos_t k = 0; k < width_at(ctx, 0); k++) {
        for (int s = 0; s < N_STATES; s++) {
            total += forward_at(&front, k, s) * back[k][s];
        }
    }

    return total;
}

/* Returns how the two passes ended.
 *
 * A total that is not finite is the backward pass overflowing: it multiplies by each row's
 * scale factor as it descends, and on a read the rates make astronomically unlikely those
 * factors are large enough that their product leaves the range of a double. That is the
 * same absence of a path a forward total of zero is, reached from the other end.
 *
 * A finite total away from one is an index error, the two passes having described
 * different sets of paths. */
static phmm_status passes_agree(const context *ctx)
{
    double total = agreement(ctx);

    if (!isfinite(total)) {
        return PHMM_NO_PATH;
    }

    return fabs(total - 1.0) < NORMALIZATION_TOLERANCE ? PHMM_OK : PHMM_UNSOUND;
}

/* Runs the backward pass. A read places at least one base, so there are always two rows
 * for the ends of it and the loop between them may be empty. */
static phmm_status backward(const context *ctx)
{
    backward_last_row(ctx);

    for (size_t i = ctx->rows - 1; i-- > 1; ) {
        backward_row(ctx, i);
    }

    backward_first_row(ctx);

    return passes_agree(ctx);
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
    if (!scratch) {
        return;
    }

    free(scratch->places);
    free(scratch->forward);
    free(scratch->terms);
    free(scratch->scale);
    free(scratch->backward);
    free(scratch->coverage);
    free(scratch->spanned);
    free(scratch->mutations);
    free(scratch);
}

/* Grows the per-row buffers. A buffer that grew is kept even where a later one
 * fails, so a scratch that fails to grow remains usable; nothing shrinks.
 *
 * Grown before the band, since the row count comes from the CIGAR while the
 * widest row is not known until the places have been written. */
static int grow_rows(phmm_scratch *scratch, size_t rows)
{
    aln_place *places;
    double    *scale;

    if (rows <= scratch->rows) {
        return 0;
    }

    places = realloc(scratch->places, rows * sizeof *places);
    scale  = realloc(scratch->scale, rows * sizeof *scale);

    if (places) {
        scratch->places = places;
    }
    if (scale) {
        scratch->scale = scale;
    }

    if (!places || !scale) {
        return -1;
    }

    scratch->rows = rows;

    return 0;
}

/* Grows the matrix buffers, which are sized by rows times widest row, so a band
 * wider than the previous read's forces a regrow however few rows this one
 * needs. */
static int grow_band(phmm_scratch *scratch, size_t rows, size_t widest)
{
    band_cell  *forward;
    band_cell  *backward;
    cell_terms *terms;

    if (rows <= scratch->matrix_rows && widest <= scratch->widest) {
        return 0;
    }

    rows   = rows   > scratch->matrix_rows ? rows   : scratch->matrix_rows;
    widest = widest > scratch->widest      ? widest : scratch->widest;

    forward  = realloc(scratch->forward, rows * widest * sizeof *forward);
    backward = realloc(scratch->backward, 2 * widest * sizeof *backward);
    terms    = realloc(scratch->terms, rows * widest * sizeof *terms);

    if (forward) {
        scratch->forward = forward;
    }
    if (backward) {
        scratch->backward = backward;
    }
    if (terms) {
        scratch->terms = terms;
    }

    if (!forward || !backward || !terms) {
        return -1;
    }

    scratch->matrix_rows = rows;
    scratch->widest      = widest;

    return 0;
}

static int grow_window(phmm_scratch *scratch, size_t window)
{
    double *coverage;
    double *spanned;
    double *mutations;

    if (window <= scratch->window) {
        return 0;
    }

    coverage  = realloc(scratch->coverage, window * sizeof *coverage);
    spanned   = realloc(scratch->spanned, window * sizeof *spanned);
    mutations = realloc(scratch->mutations, window * sizeof *mutations);

    if (coverage) {
        scratch->coverage = coverage;
    }
    if (spanned) {
        scratch->spanned = spanned;
    }
    if (mutations) {
        scratch->mutations = mutations;
    }

    if (!coverage || !spanned || !mutations) {
        return -1;
    }

    scratch->window = window;

    return 0;
}

/* ------------------------------------------------------------------------ */
/* One read                                                                  */
/* ------------------------------------------------------------------------ */

/* Returns the widest row, which is the stride every row is stored at. At least one, so
 * that an absurd band still leaves an addressable matrix; its rows are then empty and
 * the pass fails on the first. */
static hts_pos_t widest_row(const context *ctx)
{
    hts_pos_t widest = 1;

    for (size_t i = 0; i < ctx->rows; i++) {
        hts_pos_t width = width_at(ctx, i);

        if (width > widest) {
            widest = width;
        }
    }

    return widest;
}

/* Returns every position any cell can reach, plus the one before the earliest, where a
 * row's first cell lays what it pairs. Rows may widen or narrow along the read,
 * so neither end belongs to a fixed row and both are searched for. */
static extent window_of(const context *ctx)
{
    hts_pos_t first = origin_of(ctx, 0);
    hts_pos_t last  = first + width_at(ctx, 0);

    for (size_t i = 1; i < ctx->rows; i++) {
        hts_pos_t lo = origin_of(ctx, i);
        hts_pos_t hi = lo + width_at(ctx, i);

        if (lo < first) {
            first = lo;
        }
        if (hi > last) {
            last = hi;
        }
    }

    return (extent){
        .origin = first - 1,
        .len    = (size_t)(last - first) + 1,
    };
}

/* Clears the window. It is accumulated into and never assigned, so it must be cleared
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

    if (grow_rows(ctx->scratch, (size_t)read->l_qseq + 1) < 0) {
        return false;
    }

    ctx->span   = aln_places(read, ctx->scratch->places);
    ctx->rows   = (size_t)(ctx->span.end - ctx->span.begin) + 1;
    ctx->widest = widest_row(ctx);

    if (grow_band(ctx->scratch, ctx->rows, (size_t)ctx->widest) < 0) {
        return false;
    }

    ctx->window = window_of(ctx);

    if (grow_window(ctx->scratch, ctx->window.len) < 0) {
        return false;
    }

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
    phmm_status status;

    if (!prepare(&ctx)) {
        return PHMM_NO_MEMORY;
    }

    status = forward(&ctx);

    if (status != PHMM_OK) {
        return status;
    }

    status = backward(&ctx);

    if (status != PHMM_OK) {
        return status;
    }

    out->origin    = ctx.window.origin;
    out->len       = ctx.window.len;
    out->coverage  = scratch->coverage;
    out->spanned   = scratch->spanned;
    out->mutations = scratch->mutations;

    return PHMM_OK;
}
