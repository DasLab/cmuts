/* phmm.c -- a banded pair HMM.
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

/* Forward times backward sums to one on the first row within this tolerance. */
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
    double    *mismatches;
    double    *insertions;
    double    *deletions;
    double    *ends;
    size_t     rows;         /* rows places and scale are sized for */
    size_t     matrix_rows;  /* rows the forward matrix is sized for */
    size_t     widest;       /* cells a row of it holds */
    size_t     window;       /* positions the last five are sized for */
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
    const phmm_profile    *profile;
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

void phmm_build(phmm *model, const phmm_params *params)
{
    model->match_to_insertion     = 1.0 - params->extend_insertion;
    model->match_to_deletion      = 1.0 - params->extend_deletion;
    model->insertion_to_insertion = params->extend_insertion;
    model->deletion_to_deletion   = params->extend_deletion;
}

/* Returns a profile value at a 0-based base, clamped to the reference. A clamped read
 * serves a cell outside the reference, whose terms a zero emission clears, so the
 * value itself never matters. */
static double profile_at(const double *values, size_t len, hts_pos_t b)
{
    if (b < 0) {
        b = 0;
    }
    if ((size_t)b >= len) {
        b = (hts_pos_t)len - 1;
    }

    return values[b];
}

/* The three arms of the decision at the pairing of base b, which the matrix takes on
 * the edges into that pairing: continue to the next pairing, open the insertion
 * counted at b, or open the deletion counted at base b - 1. */
typedef struct {
    double match_to_match;
    double insertion_to_match;
    double deletion_to_match;
} decision;

static decision decision_at(const context *ctx, hts_pos_t b)
{
    size_t len            = ctx->ref->len;
    double open_insertion = profile_at(ctx->profile->open_insertion, len, b);
    double open_deletion  = profile_at(ctx->profile->open_deletion, len, b - 1);

    return (decision){
        .match_to_match     = 1.0 - open_insertion - open_deletion,
        .insertion_to_match = open_insertion,
        .deletion_to_match  = open_deletion,
    };
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
static double phmm_modification(double modification, bool agree, double error)
{
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
 * it that a real template difference explains, under the given modification rate. */
static cell_terms terms_from(double m, bool agree, double error)
{
    return (cell_terms){
        .emission     = agree ? agreement_chance(m, error)
                              : disagreement_chance(m, error),
        .modification = phmm_modification(m, agree, error),
    };
}

/* The parts of one row's comparisons that do not vary along the row: the row's base,
 * whether it names one, and its chance of being misread. The modification rate varies
 * per reference base, so the comparisons themselves are formed per cell. */
typedef struct {
    cell_terms neither;
    double     error;
    nuc        ours;
    bool       named;
} row_terms;

static row_terms row_terms_of(const context *ctx, size_t i)
{
    int32_t query = ctx->span.begin + (int32_t)i - 1;
    nuc     ours  = nuc_from_read(ctx->read->seq, query);

    return (row_terms){
        .neither = { .emission = UNINFORMATIVE, .modification = 0.0 },
        .error   = error_at(ctx, query),
        .ours    = ours,
        .named   = nuc_is_base(ours),
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

    if (!nuc_is_base(theirs) || !row->named) {
        return row->neither;
    }

    return terms_from(ctx->profile->modification[j - 1], theirs == row->ours,
                      row->error);
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

/* The row above one forward row and the transition weights out of it that hold along
 * the row, with the row above's scale factor and the insertion emission folded in. The
 * decision weights vary per cell and are taken at each cell instead. The folded scale
 * factor is why this is the one place a forward row is read without forward_at. */
typedef struct {
    const band_cell *above;
    hts_pos_t        width;
    hts_pos_t        shift;
    double           scale;
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
        .scale                  = scale,
        .match_to_insertion     = model->match_to_insertion * scale
                                * UNINFORMATIVE,
        .insertion_to_insertion = model->insertion_to_insertion * scale
                                * UNINFORMATIVE,
    };
}

/* Returns cell k's match state, stepped from the row above through the decision at
 * the cell's own base. */
static double paired_from(const descent *step, const decision *into, hts_pos_t k,
                          double emission)
{
    hts_pos_t diagonal = k - 1 + step->shift;  /* a position back, one row up */

    if (!within(diagonal, step->width)) {
        return 0.0;
    }

    return (into->match_to_match * step->scale
              * step->above[diagonal][STATE_MATCH]
          + into->insertion_to_match * step->scale
              * step->above[diagonal][STATE_INSERTION]
          + into->deletion_to_match * step->scale
              * step->above[diagonal][STATE_DELETION])
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

/* Returns whether row i can hold a deletion. */
static bool deletions_live(const context *ctx, size_t i)
{
    return i > 0 && i + 1 < ctx->rows;
}

/* Returns whether row i can hold an insertion. */
static bool insertions_live(const context *ctx, size_t i)
{
    return i > 1 && i + 1 < ctx->rows;
}

/* Fills row i of the forward pass and returns its unscaled total. Valid only
 * for a row with a row above it. */
static double forward_row(const context *ctx, size_t i)
{
    const phmm *model      = ctx->model;
    descent     step       = descent_into(ctx, i);
    band_cell  *row        = row_of(ctx, i);
    cell_terms *terms      = terms_of(ctx, i);
    row_terms   each       = row_terms_of(ctx, i);
    hts_pos_t   width      = width_at(ctx, i);
    bool        deletions  = deletions_live(ctx, i);
    bool        insertions = insertions_live(ctx, i);
    /* The cell to the left, held in locals so each step of the deletion chain
     * does not wait on the preceding store. */
    double      left_match    = 0.0;
    double      left_deletion = 0.0;
    /* One running sum per state keeps the addition chains short. */
    double      total_paired   = 0.0;
    double      total_inserted = 0.0;
    double      total_deleted  = 0.0;

    for (hts_pos_t k = 0; k < width; k++) {
        hts_pos_t j = position_of(ctx, i, k);
        decision  into = decision_at(ctx, j - 1);  /* the base the cell pairs */
        double    paired, inserted, deleted;

        terms[k] = terms_at(ctx, &each, j);

        paired   = paired_from(&step, &into, k, terms[k].emission);
        inserted = insertions ? inserted_from(&step, k) : 0.0;
        deleted  = deletions && k > 0
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

/* The window fields, each advanced to where the row's first cell enters
 * the window, so a cell addresses its positions by its own index with no
 * bounds check. window_of guarantees that every position a row can address
 * lies inside the window. The band is not clamped, so a row near either end
 * of the reference addresses positions outside it; those take only zero
 * contributions. */
typedef struct {
    double *coverage;
    double *mismatches;
    double *insertions;
    double *deletions;
    double *ends;
} landing;

static landing landing_of(const context *ctx, size_t i)
{
    phmm_scratch *scratch = ctx->scratch;
    size_t        at      = (size_t)(origin_of(ctx, i) - 1
                                   - ctx->window.origin);

    return (landing){
        .coverage   = scratch->coverage + at,
        .mismatches = scratch->mismatches + at,
        .insertions = scratch->insertions + at,
        .deletions  = scratch->deletions + at,
        .ends       = scratch->ends + at,
    };
}

/* The state of one row's accumulation into the window. Cells are accumulated
 * right to left, the order the backward pass forms them in. The pairing of a
 * position and the deletion of it are held by the cell to its right, so a
 * position is complete only once the cell to its left is reached. */
typedef struct {
    scaled_row        front;
    const cell_terms *terms;
    landing           at;
    /* Where the row's pairings are also the read's 5'-most, which no decision
     * follows; every other row leaves this NULL. */
    double           *ends;
    /* The pending contribution to the position the next cell completes. */
    double            coverage;
    double            mismatches;
    double            deletions;
} accumulation;

static accumulation accumulation_of(const context *ctx, size_t i)
{
    landing at = landing_of(ctx, i);

    return (accumulation){
        .front = scaled_row_of(ctx, i),
        .terms = terms_of(ctx, i),
        .at    = at,
        /* Row 1 pairs the read's 5'-most placed base. */
        .ends  = i == 1 ? at.ends : NULL,
    };
}

/* Accumulates cell k of the row into the window. The two runs a cell can open begin
 * on the decision at the pairing below, whose weights the caller passes in. */
static void accumulate_cell(accumulation *acc, hts_pos_t k, const decision *below,
                            const double *back, double pairing)
{
    double matched  = forward_at(&acc->front, k, STATE_MATCH);
    double skipped  = forward_at(&acc->front, k, STATE_DELETION);
    double carried  = forward_at(&acc->front, k, STATE_INSERTION);
    double paired   = matched * back[STATE_MATCH];
    double deleted  = below->deletion_to_match * skipped * pairing;
    double inserted = below->insertion_to_match * carried * pairing;

    /* Accumulated from the previous call, except the insertion. */
    acc->at.coverage[k + 1]   += acc->coverage;
    acc->at.mismatches[k + 1] += acc->mismatches;
    acc->at.deletions[k + 1]  += acc->deletions;
    acc->at.insertions[k + 1] += inserted;

    /* Read by the next call. */
    acc->coverage   = paired;
    acc->mismatches = paired * acc->terms[k].modification;
    acc->deletions  = deleted;

    /* Written at the index the pending values of this cell land at. */
    if (acc->ends) {
        acc->ends[k] += paired;
    }
}

/* Writes the pending contribution after the leftmost cell, which completes
 * the row's last open position. */
static void accumulate_end(const accumulation *acc)
{
    acc->at.coverage[0]   += acc->coverage;
    acc->at.mismatches[0] += acc->mismatches;
    acc->at.deletions[0]  += acc->deletions;
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

/* Forms a cell's three states from the transitions out of it. The steps into the
 * pairing below carry the decision at that pairing's base. */
static void backward_cell(const phmm *model, const decision *below, bool insertions,
                          double pairing, double inserted, double deleted,
                          double *cell)
{
    cell[STATE_MATCH] = below->match_to_match     * pairing
                      + model->match_to_insertion * inserted
                      + model->match_to_deletion  * deleted;

    cell[STATE_INSERTION] = insertions
                          ? below->insertion_to_match     * pairing
                          + model->insertion_to_insertion * inserted
                          : 0.0;

    cell[STATE_DELETION] = below->deletion_to_match    * pairing
                         + model->deletion_to_deletion * deleted;
}

/* Fills the last row of the backward pass and accumulates it. The alignment ends on
 * this row, and it ends on a pairing, so a match finishes with chance one and neither
 * run can be open. */
static void backward_last_row(const context *ctx)
{
    size_t       i   = ctx->rows - 1;
    band_cell   *row = backward_row_of(ctx, i);
    accumulation acc = accumulation_of(ctx, i);
    double       cell[N_STATES] = {
        [STATE_MATCH]     = 1.0,
        [STATE_INSERTION] = 0.0,
        [STATE_DELETION]  = 0.0,
    };

    for (hts_pos_t k = width_at(ctx, i); k-- > 0; ) {
        decision below = decision_at(ctx, position_of(ctx, i, k));

        row[k][STATE_MATCH]     = cell[STATE_MATCH];
        row[k][STATE_INSERTION] = cell[STATE_INSERTION];
        row[k][STATE_DELETION]  = cell[STATE_DELETION];

        accumulate_cell(&acc, k, &below, cell, 0.0);
    }

    accumulate_end(&acc);
}

/* Fills row i of the backward pass and accumulates it. Valid only for a row
 * with a row above and a row below. */
static void backward_row(const context *ctx, size_t i)
{
    const phmm  *model      = ctx->model;
    band_cell   *row        = backward_row_of(ctx, i);
    ascent       below      = ascent_into(ctx, i);
    accumulation acc        = accumulation_of(ctx, i);
    hts_pos_t    width      = width_at(ctx, i);
    bool         insertions = insertions_live(ctx, i);
    /* The cell to the right, held in a local so each step of the deletion
     * chain does not wait on the preceding store. */
    double       right_deletion = 0.0;

    for (hts_pos_t k = width; k-- > 0; ) {
        decision into     = decision_at(ctx, position_of(ctx, i, k));
        double   pairing  = pairing_below(&below, k);
        double   inserted = inserted_below(&below, k);
        double   deleted  = k + 1 < width ? right_deletion : 0.0;
        double   cell[N_STATES];

        backward_cell(model, &into, insertions, pairing, inserted, deleted, cell);

        row[k][STATE_MATCH]     = cell[STATE_MATCH];
        row[k][STATE_INSERTION] = cell[STATE_INSERTION];
        row[k][STATE_DELETION]  = cell[STATE_DELETION];

        accumulate_cell(&acc, k, &into, cell, pairing);

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
        decision into     = decision_at(ctx, position_of(ctx, 0, k));
        double   pairing  = pairing_below(&below, k);
        double   inserted = inserted_below(&below, k);

        row[k][STATE_MATCH] = into.match_to_match       * pairing
                            + model->match_to_insertion * inserted;

        /* Neither run reaches the row before the first placed base. */
        row[k][STATE_INSERTION] = 0.0;
        row[k][STATE_DELETION]  = 0.0;
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
    free(scratch->mismatches);
    free(scratch->insertions);
    free(scratch->deletions);
    free(scratch->ends);
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
    double *mismatches;
    double *insertions;
    double *deletions;
    double *ends;

    if (window <= scratch->window) {
        return 0;
    }

    coverage   = realloc(scratch->coverage, window * sizeof *coverage);
    mismatches = realloc(scratch->mismatches, window * sizeof *mismatches);
    insertions = realloc(scratch->insertions, window * sizeof *insertions);
    deletions  = realloc(scratch->deletions, window * sizeof *deletions);
    ends       = realloc(scratch->ends, window * sizeof *ends);

    if (coverage) {
        scratch->coverage = coverage;
    }
    if (mismatches) {
        scratch->mismatches = mismatches;
    }
    if (insertions) {
        scratch->insertions = insertions;
    }
    if (deletions) {
        scratch->deletions = deletions;
    }
    if (ends) {
        scratch->ends = ends;
    }

    if (!coverage || !mismatches || !insertions || !deletions || !ends) {
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
    memset(scratch->mismatches, 0, len * sizeof *scratch->mismatches);
    memset(scratch->insertions, 0, len * sizeof *scratch->insertions);
    memset(scratch->deletions, 0, len * sizeof *scratch->deletions);
    memset(scratch->ends, 0, len * sizeof *scratch->ends);
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

phmm_status phmm_run(const phmm *model, const phmm_profile *profile,
                     const phred *quality, const cm_bam_record *read,
                     const cm_fasta_record *ref, const int *half,
                     phmm_scratch *scratch, phmm_window *out)
{
    context ctx = {
        .model   = model,
        .profile = profile,
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

    out->origin     = ctx.window.origin;
    out->len        = ctx.window.len;
    out->coverage   = scratch->coverage;
    out->mismatches = scratch->mismatches;
    out->insertions = scratch->insertions;
    out->deletions  = scratch->deletions;
    out->ends       = scratch->ends;

    return PHMM_OK;
}
