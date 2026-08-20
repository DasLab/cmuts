/* phmm.c -- a banded pair HMM over one read and the reference it was placed on.
 *
 * Three states: M pairs a read base with a reference base, I consumes a read
 * base with no reference base, D consumes a reference base with no read base.
 * Both indices are prefix lengths, so cell (i, j) is the first i placed read
 * bases against the first j reference bases, with the state giving the last
 * operation. A cell in M or I compares reference base j - 1 with read base
 * i - 1.
 *
 * Each read base carries a band of reference positions: what the CIGAR path
 * crosses on that row, plus a half-width on each side. Every probability is
 * conditional on the alignment staying inside the band. The band marginalizes
 * local ambiguity but does not correct a misplaced read. A band of 0
 * marginalizes over the CIGAR path only.
 *
 * Both passes divide each row by the forward row's total. A posterior is the
 * product of the two passes, with no separate normalizer and no logarithms.
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

/* Emission of a comparison that carries no information: an inserted base, or a
 * base one side does not hold. */
#define UNINFORMATIVE (1.0 / NUC_BASES)

/* A misread or modified base is equally likely to be any base the reference
 * does not hold. */
#define OTHER_BASES ((double)(NUC_BASES - 1))

/* Forward times backward sums to one on the first row. A finite departure past
 * this tolerance indicates an index error, not rounding, and the run stops. */
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

/* The inputs and derived layout for one read.
 *
 * The caller's inputs come first and are fixed for the read. prepare() derives
 * the rest and is the only writer, so every other function takes a const
 * context. */
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

/* Returns the chance a read base agrees with the base it was templated from.
 * A base can disagree because the template was modified or because the base
 * was misread. Agreement plus three disagreements sums to one. */
static double agreement_chance(double modification, double error)
{
    return (1.0 - modification) * (1.0 - error)
         + modification * error / OTHER_BASES;
}

/* Returns the chance a read base differs from the base it was templated
 * from. */
static double disagreement_chance(double modification, double error)
{
    return (1.0 - modification) * error / OTHER_BASES
         + modification * (1.0 - error / OTHER_BASES) / OTHER_BASES;
}

/* Returns the chance the template differed here, given the base as read. An
 * agreeing base may have been modified and misread back into agreement. A
 * disagreeing base may be an unmodified base misread, so a poorly read
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

/* Returns the chance the base at the given query offset was misread, or zero
 * when the record stores no qualities. */
static double error_at(const context *ctx, int32_t query)
{
    return ctx->read->qual
         ? phred_error(ctx->quality, ctx->read->qual[query])
         : 0.0;
}

/* Returns the emission of an agreeing or disagreeing pairing, and the part of
 * it that a real template difference explains. */
static cell_terms terms_from(const context *ctx, bool agree, double error)
{
    double m = ctx->model->params.modification;

    return (cell_terms){
        .emission     = agree ? agreement_chance(m, error)
                              : disagreement_chance(m, error),
        .modification = phmm_modification(ctx->model, agree, error),
    };
}

/* The three comparisons the cells of one row can hold, precomputed on entry
 * to the row: the reference base agrees with the row's base, differs from it,
 * or is missing. Where the read holds no base the three are equal. */
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

/* Returns the comparison for the cell at reference prefix length j, which
 * pairs the row's base with reference base j - 1. A position past either end
 * of the reference has no base, so its emission is zero. An ambiguous base is
 * a real base of unknown identity and takes the uninformative emission. */
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

/* ------------------------------------------------------------------------ */
/* The band                                                                  */
/* ------------------------------------------------------------------------ */

/* Returns a forward row as stored, which is not the value it represents. See
 * scaled_row. */
static band_cell *row_of(const context *ctx, size_t i)
{
    return ctx->scratch->forward + i * (size_t)ctx->widest;
}

/* Returns the same row, const. A pointer to an array takes its own qualifier
 * only from C23, so the cast is here and in read_backward_row_of, not at every
 * call site. */
static const band_cell *read_row_of(const context *ctx, size_t i)
{
    return (const band_cell *)row_of(ctx, i);
}

/* A forward row and the factor that restores its true values.
 *
 * The forward pass stores each row divided by the totals of the rows above
 * it, which keeps the values from underflowing over long reads. A stored cell
 * is not a probability on its own. Read cells through forward_at, which
 * applies the factor. */
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

/* Returns the count of reference bases the CIGAR deletes directly after the
 * row's base. */
static hts_pos_t skip_at(const context *ctx, size_t i)
{
    return ctx->scratch->places[i].last - ctx->scratch->places[i].first;
}

/* Returns the cell count of row i: the stretch the CIGAR path crosses on it,
 * plus a half-width on each side. */
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

/* Returns the offset that maps a cell index of one row to the same position on
 * another row, which is the distance between their first cells. */
static hts_pos_t shift_between(const context *ctx, size_t from, size_t to)
{
    return origin_of(ctx, to) - origin_of(ctx, from);
}

/* Returns whether k is a cell index of a row with the given width. */
static bool within(hts_pos_t k, hts_pos_t width)
{
    return k >= 0 && k < width;
}

/* ------------------------------------------------------------------------ */
/* Forward                                                                   */
/* ------------------------------------------------------------------------ */

/* Fills the first row of the forward pass and returns its total. The
 * alignment starts at each position of the row with equal chance, and never
 * in a deletion. */
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

/* The row above one forward row and the transition weights out of it, with
 * the row above's scale factor and the insertion emission folded into the
 * weights. The folded scale factor is why this is the one place a forward row
 * is read without forward_at. */
typedef struct {
    const band_cell *above;
    hts_pos_t        width;
    hts_pos_t        shift;
    double           match_to_match;
    double           insertion_to_match;
    double           deletion_to_match;
    double           match_to_insertion;
    double           insertion_to_insertion;
} descent;

static descent descent_into(const context *ctx, size_t i)
{
    const phmm *model = ctx->model;
    double      scale = ctx->scratch->scale[i - 1];

    return (descent){
        .above                  = read_row_of(ctx, i - 1),
        .width                  = width_at(ctx, i - 1),
        .shift                  = shift_between(ctx, i - 1, i),
        .match_to_match         = model->match_to_match * scale,
        .insertion_to_match     = model->insertion_to_match * scale,
        .deletion_to_match      = model->deletion_to_match * scale,
        .match_to_insertion     = model->match_to_insertion * scale
                                * UNINFORMATIVE,
        .insertion_to_insertion = model->insertion_to_insertion * scale
                                * UNINFORMATIVE,
    };
}

/* Returns cell k's match state, stepped from the row above. */
static double paired_from(const descent *step, hts_pos_t k, double emission)
{
    hts_pos_t diagonal = k - 1 + step->shift;  /* a position back, one row up */

    if (!within(diagonal, step->width)) {
        return 0.0;
    }

    return (step->match_to_match     * step->above[diagonal][STATE_MATCH]
          + step->insertion_to_match * step->above[diagonal][STATE_INSERTION]
          + step->deletion_to_match  * step->above[diagonal][STATE_DELETION])
         * emission;
}

/* Returns cell k's insertion state, stepped from the row above. */
static double inserted_from(const descent *step, hts_pos_t k)
{
    hts_pos_t straight = k + step->shift;      /* this position, one row up */

    if (!within(straight, step->width)) {
        return 0.0;
    }

    return step->match_to_insertion     * step->above[straight][STATE_MATCH]
         + step->insertion_to_insertion * step->above[straight][STATE_INSERTION];
}

/* Returns a cell's deletion state, stepped from the cell to its left. */
static double deleted_from(const phmm *model, double left_match,
                           double left_deletion)
{
    return model->match_to_deletion    * left_match
         + model->deletion_to_deletion * left_deletion;
}

/* Returns whether row i can hold a deletion. The first and last rows cannot:
 * the alignment neither starts nor ends in one. */
static bool deletions_live(const context *ctx, size_t i)
{
    return i > 0 && i + 1 < ctx->rows;
}

/* Fills row i of the forward pass and returns its unscaled total. Valid only
 * for a row with a row above it. */
static double forward_row(const context *ctx, size_t i)
{
    const phmm *model = ctx->model;
    descent     step  = descent_into(ctx, i);
    band_cell  *row   = row_of(ctx, i);
    cell_terms *terms = terms_of(ctx, i);
    row_terms   each  = row_terms_of(ctx, i);
    hts_pos_t   width = width_at(ctx, i);
    bool        live  = deletions_live(ctx, i);
    /* The cell to the left, held in locals so each step of the deletion chain
     * does not wait on the preceding store. */
    double      left_match    = 0.0;
    double      left_deletion = 0.0;
    /* One running sum per state keeps the addition chains short. */
    double      total_paired   = 0.0;
    double      total_inserted = 0.0;
    double      total_deleted  = 0.0;

    for (hts_pos_t k = 0; k < width; k++) {
        double paired, inserted, deleted;

        terms[k] = terms_at(ctx, &each, position_of(ctx, i, k));

        paired   = paired_from(&step, k, terms[k].emission);
        inserted = inserted_from(&step, k);
        deleted  = live && k > 0
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

/* Checks row i's total and stores its reciprocal as the row's scale factor.
 * A zero or subnormal total means the band admits no path for this read,
 * which is not a fault. A negative or non-finite total is a fault. */
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

/* Runs the forward pass. The first row is filled directly, since it has no row
 * above to step from. */
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

/* Each reference position collects three totals: coverage, span, and
 * mutations. A pairing covers and spans the position it pairs, and adds as
 * mutation the part of its posterior that a template modification explains.
 * A deletion spans every position it passes over, adds no coverage, and counts as
 * one mutation at the end of its run, since reverse transcription reads the
 * template from the 3' end. An insertion counts as a mutation at the position
 * it precedes and adds no coverage; it spans the same weighted amount it counts
 * as mutation, so a weight of zero removes it entirely. */

/* The three window fields, each advanced to where the row's first cell enters
 * the window, so a cell addresses its positions by its own index with no
 * bounds check. window_of guarantees that every position a row can address
 * lies inside the window. The band is not clamped, so a row near either end
 * of the reference addresses positions outside it; those take only zero
 * contributions. */
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

/* The weight of one event of each mutation kind, with every factor of the
 * event's posterior that is constant along the row folded in. Valid only for
 * a row with a row above it. */
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

/* The state of one row's accumulation into the window. Cells are accumulated
 * right to left, the order the backward pass forms them in. A position is
 * complete only when the cell to its left is reached, since that cell can
 * carry an insertion into it. Valid only for a row with a row above it. */
typedef struct {
    scaled_row        front;
    scaled_row        above;
    const cell_terms *terms;
    weighing          weight;
    landing           at;
    hts_pos_t         up;           /* to the same position, one row up */
    hts_pos_t         above_width;
    /* The pending contribution to the position the next cell completes. */
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
        .up          = shift_between(ctx, i - 1, i),
        .above_width = width_at(ctx, i - 1),
    };
}

/* Accumulates cell k of the row: writes the position the cell completes, and
 * holds the cell's own contribution for the position to its left. back is the
 * cell's backward states. pairing is the backward pairing of the row below,
 * which a deletion run closes against. */
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

    acc->coverage  = paired;
    acc->spanned   = paired + passed;
    acc->mutations = acc->weight.substitution * paired
                   * acc->terms[k].modification
                   + acc->weight.deletion * skipped * pairing;
}

/* Writes the pending contribution after the leftmost cell, which completes
 * the row's last open position. */
static void accumulate_end(const accumulation *acc)
{
    acc->at.coverage[0]  += acc->coverage;
    acc->at.spanned[0]   += acc->spanned;
    acc->at.mutations[0] += acc->mutations;
}

/* ------------------------------------------------------------------------ */
/* Backward                                                                  */
/* ------------------------------------------------------------------------ */

/* A backward cell holds the chance of every way the alignment can finish
 * from that cell. Each step into a row multiplies in that row's forward scale
 * factor, so a forward cell times a backward cell is a posterior. Each row is
 * accumulated into the window as it is formed. */

/* The row below one backward row: its cells, comparisons, scale factor,
 * width, and the shift between the rows. */
typedef struct {
    const band_cell  *cell;
    const cell_terms *terms;
    double            scale;
    hts_pos_t         width;
    hts_pos_t         shift;
} ascent;

static ascent ascent_into(const context *ctx, size_t i)
{
    return (ascent){
        .cell  = read_backward_row_of(ctx, i + 1),
        .terms = terms_of(ctx, i + 1),
        .scale = ctx->scratch->scale[i + 1],
        .width = width_at(ctx, i + 1),
        .shift = shift_between(ctx, i, i + 1),
    };
}

/* Returns the chance of finishing through a pairing of the next read base,
 * one position on and one row down from cell k. */
static double pairing_below(const ascent *below, hts_pos_t k)
{
    hts_pos_t diagonal = k + 1 - below->shift;

    if (!within(diagonal, below->width)) {
        return 0.0;
    }

    return below->terms[diagonal].emission
         * below->cell[diagonal][STATE_MATCH] * below->scale;
}

/* Returns the chance of finishing through an insertion of the next read base,
 * at the same position one row down from cell k. */
static double inserted_below(const ascent *below, hts_pos_t k)
{
    hts_pos_t straight = k - below->shift;

    if (!within(straight, below->width)) {
        return 0.0;
    }

    return UNINFORMATIVE
         * below->cell[straight][STATE_INSERTION] * below->scale;
}

/* Forms a cell's three states from the transitions out of it. */
static void backward_cell(const phmm *model, double pairing, double inserted,
                          double deleted, double *cell)
{
    cell[STATE_MATCH] = model->match_to_match     * pairing
                      + model->match_to_insertion * inserted
                      + model->match_to_deletion  * deleted;

    cell[STATE_INSERTION] = model->insertion_to_match     * pairing
                          + model->insertion_to_insertion * inserted;

    cell[STATE_DELETION] = model->deletion_to_match    * pairing
                         + model->deletion_to_deletion * deleted;
}

/* Fills the last row of the backward pass and accumulates it. The alignment
 * ends on this row, so a match or insertion finishes with chance one and a
 * deletion cannot occur. */
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

/* Fills row i of the backward pass and accumulates it. Valid only for a row
 * with a row above and a row below. */
static void backward_row(const context *ctx, size_t i)
{
    const phmm  *model = ctx->model;
    band_cell   *row   = backward_row_of(ctx, i);
    ascent       below = ascent_into(ctx, i);
    accumulation acc   = accumulation_of(ctx, i);
    hts_pos_t    width = width_at(ctx, i);
    /* The cell to the right, held in a local so each step of the deletion
     * chain does not wait on the preceding store. */
    double       right_deletion = 0.0;

    for (hts_pos_t k = width; k-- > 0; ) {
        double pairing  = pairing_below(&below, k);
        double inserted = inserted_below(&below, k);
        double deleted  = k + 1 < width ? right_deletion : 0.0;
        double cell[N_STATES];

        backward_cell(model, pairing, inserted, deleted, cell);

        row[k][STATE_MATCH]     = cell[STATE_MATCH];
        row[k][STATE_INSERTION] = cell[STATE_INSERTION];
        row[k][STATE_DELETION]  = cell[STATE_DELETION];

        accumulate_cell(&acc, k, cell, pairing);

        right_deletion = cell[STATE_DELETION];
    }

    accumulate_end(&acc);
}

/* Fills the first row of the backward pass. The row precedes the first placed
 * base, so it adds no value to the window and is written only for
 * passes_agree to check. */
static void backward_first_row(const context *ctx)
{
    const phmm *model = ctx->model;
    band_cell  *row   = backward_row_of(ctx, 0);
    ascent      below = ascent_into(ctx, 0);
    hts_pos_t   width = width_at(ctx, 0);

    for (hts_pos_t k = 0; k < width; k++) {
        double pairing  = pairing_below(&below, k);
        double inserted = inserted_below(&below, k);

        row[k][STATE_MATCH] = model->match_to_match     * pairing
                            + model->match_to_insertion * inserted;

        row[k][STATE_INSERTION] = model->insertion_to_match     * pairing
                                + model->insertion_to_insertion * inserted;

        /* No deletion occurs before the first placed base. */
        row[k][STATE_DELETION] = 0.0;
    }
}

/* Returns the first row's total of forward times backward, which is one when
 * the two passes describe the same set of paths. The first row carries the
 * whole backward chain and every forward scale factor, so this checks both
 * passes end to end. */
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

/* Maps the agreement total to a status. A total that is not finite means the
 * backward pass overflowed on the scale factors, which happens only when the
 * band admits no plausible path. A finite total away from one is an index
 * error. */
static phmm_status passes_agree(const context *ctx)
{
    double total = agreement(ctx);

    if (!isfinite(total)) {
        return PHMM_NO_PATH;
    }

    return fabs(total - 1.0) < NORMALIZATION_TOLERANCE ? PHMM_OK : PHMM_UNSOUND;
}

/* Runs the backward pass. A read places at least one base, so there are always
 * two rows for the ends, and the loop between them can be empty. */
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

/* Grows the per-row buffers. A buffer that grew is kept even when a later one
 * fails, so a failed grow leaves the scratch usable; the buffers never shrink. This
 * runs before grow_band, since the row count comes from the CIGAR while the
 * widest row is known only after the places are written. */
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

/* Grows the matrix buffers, which are sized as rows times the widest row. A
 * band wider than the previous read's forces a regrow even when this read
 * needs fewer rows. */
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

/* Returns the width of the widest row, and at least one, which is the stride
 * every row is stored at. */
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

/* Returns the smallest reference range that holds every position any row can
 * write. */
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

/* Clears the window. It is accumulated into and never assigned, so it must be
 * cleared between reads. */
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

void phmm_window_bounds(const phmm_window *window, size_t len,
                        size_t *from, size_t *to)
{
    hts_pos_t low  = window->origin < 0 ? -window->origin : 0;
    hts_pos_t high = (hts_pos_t)len - window->origin;

    *from = (size_t)low;
    *to   = high < 0 ? 0 : (size_t)high;

    if (*to > window->len) {
        *to = window->len;
    }
    if (*from > *to) {
        *from = *to;
    }
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
