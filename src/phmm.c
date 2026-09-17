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

/* The three arms of the decision at the pairing of base b, which the matrix takes on
 * the edges into that pairing: continue to the next pairing, open the insertion
 * counted at b, or open the deletion counted at base b - 1. */
typedef struct {
    double match_to_match;
    double insertion_to_match;
    double deletion_to_match;
} decision;

/* What a column's base is to every row: one of the named bases, an ambiguity code, or
 * absent past either end of the reference. The values below NUC_COUNT are the bases as
 * nuc names them, so that a code indexes a table of comparisons directly. */
#define COLUMN_ABSENT NUC_COUNT
#define COLUMN_KINDS  (NUC_COUNT + 1)

/* The decision at the base's pairing, its modification rate, and what its base is. An
 * absent base takes the rate of the nearest base, so that a row crossing the end of the
 * reference sees one rate throughout. */
struct phmm_column {
    decision decision;
    double   modification;
    uint8_t  code;
};

/* Where one row of the band lies: the reference prefix length its first cell stands
 * for, and its cell count. */
typedef struct {
    hts_pos_t origin;
    hts_pos_t width;
} row_band;

/* The comparisons of a read base of one quality against a reference base it agrees
 * with and one it differs from, under one modification rate. Kept across reads: the
 * qualities a read carries recur, and the rate is one value along most of a reference,
 * so the divisions that form the comparisons run once for each quality rather than once
 * for each row. */
typedef struct {
    cell_terms agree;
    cell_terms disagree;
    double     modification;  /* the rate the two hold for; NaN before any is kept,
                                 which no rate compares equal to */
} quality_terms;

/* A slot for each quality a byte can hold, and one for a read storing no qualities. */
#define QUALITY_SLOTS (PHRED_MAX + 2)
#define NO_QUALITY    (PHRED_MAX + 1)

/* Rows are stored at the widest row's stride, so a row is located by
 * multiplication. Only the loops are ragged. */
struct phmm_scratch {
    aln_place *places;     /* one per placed read base, and one before them */
    row_band  *bands;      /* one per row */
    double    *factor;     /* what each forward row is stored multiplied by */
    double     restore;    /* the reciprocal of the last forward row's total */
    band_cell *forward;    /* rows, each of widest cells */
    cell_terms *terms;     /* one per cell of it */
    /* Two rows suffice: only the current row and the one below it are read. */
    band_cell *backward;
    phmm_position *window; /* the values returned to the caller, one per position */
    quality_terms by_quality[QUALITY_SLOTS];  /* the comparisons of each quality */
    size_t     rows;         /* rows places, bands and factor are sized for */
    size_t     matrix_rows;  /* rows the forward matrix is sized for */
    size_t     widest;       /* cells a row of it holds */
    size_t     window_len;   /* positions window is sized for */
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
    const phmm_model      *model;
    const phmm_rates      *rates;   /* the model's */
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

void phmm_rates_set_transitions(phmm_rates *rates, const phmm_uniform_rates *uniform)
{
    rates->match_to_insertion     = 1.0 - uniform->extend_insertion;
    rates->match_to_deletion      = 1.0 - uniform->extend_deletion;
    rates->insertion_to_insertion = uniform->extend_insertion;
    rates->deletion_to_deletion   = uniform->extend_deletion;
}

/* Returns a 0-based base clamped to a reference of len bases. */
static hts_pos_t clamp_base(size_t len, hts_pos_t b)
{
    if (b < 0) {
        b = 0;
    }
    if ((size_t)b >= len) {
        b = (hts_pos_t)len - 1;
    }

    return b;
}

/* Returns a per-base rate at a 0-based base, clamped to the reference. A clamped read
 * serves a cell outside the reference, whose terms a zero emission clears, so the
 * value itself never matters. */
static double rate_at(const double *values, size_t len, hts_pos_t b)
{
    return values[clamp_base(len, b)];
}

static decision decision_at(const phmm_rates *rates, size_t len, hts_pos_t b)
{
    double open_insertion = rate_at(rates->open_insertion, len, b);
    double open_deletion  = rate_at(rates->open_deletion, len, b - 1);

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

/* Returns the quality slot of the base at the given query offset. */
static int quality_at(const context *ctx, int32_t query)
{
    return ctx->read->qual ? ctx->read->qual[query] : NO_QUALITY;
}

/* Returns the chance a base of the given quality slot was misread, which is zero for a
 * read storing no qualities. */
static double error_of(const context *ctx, int quality)
{
    return quality == NO_QUALITY ? 0.0
                                 : phred_error(ctx->quality, (uint8_t)quality);
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

/* Forms the comparisons of one quality for one modification rate. */
static void quality_terms_keep(quality_terms *terms, double modification, double error)
{
    terms->agree        = terms_from(modification, true, error);
    terms->disagree     = terms_from(modification, false, error);
    terms->modification = modification;
}

/* The comparisons of one row's base against each kind of column. An absent base has no
 * comparison, so its emission is zero. An ambiguous base on either side is a real base
 * of unknown identity and takes the uninformative emission. The comparisons against
 * named bases depend on the modification rate, which varies per reference base, so
 * they are drawn from the row's quality slot for one rate at a time and kept while the
 * cells along the row read that rate. */
typedef struct {
    cell_terms     by_code[COLUMN_KINDS];
    double         kept;      /* the modification rate the named comparisons hold for;
                                 NaN before any is kept */
    double         error;
    quality_terms *quality;   /* the slot of the row's base quality */
    nuc            ours;
} row_terms;

static row_terms row_terms_of(const context *ctx, size_t i)
{
    int32_t   query   = ctx->span.begin + (int32_t)i - 1;
    int       quality = quality_at(ctx, query);
    row_terms row     = {
        .kept    = (double)NAN,
        .error   = error_of(ctx, quality),
        .quality = &ctx->scratch->by_quality[quality],
        .ours    = nuc_from_read(ctx->read->seq, query),
    };

    row.by_code[COLUMN_ABSENT] = (cell_terms){ .emission = 0.0, .modification = 0.0 };
    row.by_code[NUC_N]         = (cell_terms){ .emission = UNINFORMATIVE,
                                               .modification = 0.0 };

    return row;
}

/* Forms the row's comparisons against the named bases for one modification rate,
 * drawing on the row's quality slot and refilling it where it holds another rate. */
static void row_terms_keep(row_terms *row, double modification)
{
    quality_terms *quality = row->quality;

    if (modification != quality->modification) {
        quality_terms_keep(quality, modification, row->error);
    }

    for (nuc theirs = NUC_A; theirs < NUC_COUNT; theirs++) {
        row->by_code[theirs] = !nuc_is_base(row->ours) ? row->by_code[NUC_N]
                             : theirs == row->ours     ? quality->agree
                                                       : quality->disagree;
    }

    row->kept = modification;
}

/* Returns the comparison of the row's base with the base of one column. */
static cell_terms terms_at(row_terms *row, const phmm_column *col)
{
    if (col->modification != row->kept) {
        row_terms_keep(row, col->modification);
    }

    return row->by_code[col->code];
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

/* The scaling of the passes.
 *
 * The forward pass stores each row multiplied by a factor of its own, which keeps the
 * values from underflowing over long reads, so a stored cell is not a probability on
 * its own. The factor is formed from the rows two and more above, not from the row
 * directly above: taking that row's total would put its summation and division on the
 * chain every cell of the row waits on, whereas the row two above is complete a whole
 * row earlier. The factor of row i is the reciprocal of the total of row i - 2 times
 * the factor of row i - 1, and one for the first two rows. It undoes both, so the total
 * of row i comes out near the product of the decay of rows i - 1 and i, whatever the
 * rows before them did. A factor taken from the total of row i - 2 alone would leave
 * each total the one above divided by the one before, which oscillates and grows.
 *
 * The backward pass folds the same factor into each step into a row. The product of a
 * stored forward cell and a stored backward cell is then the posterior divided by the
 * total of the last row, whose reciprocal forward_at applies. */

/* Returns the factor row i was stored with. */
static double factor_into(const context *ctx, size_t i)
{
    return ctx->scratch->factor[i];
}

/* A forward row and the factor that restores its true values. Read cells through
 * forward_at, which applies the factor. */
typedef struct {
    const band_cell *cell;
    double           scale;
} scaled_row;

static scaled_row scaled_row_of(const context *ctx, size_t i)
{
    return (scaled_row){
        .cell  = read_row_of(ctx, i),
        .scale = ctx->scratch->restore,
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

/* Returns where row i lies, from the CIGAR and the band: the stretch the CIGAR path
 * crosses on it, plus a half-width on each side. */
static row_band band_of(const context *ctx, size_t i)
{
    return (row_band){
        .origin = ctx->scratch->places[i].first - ctx->half[i],
        .width  = skip_at(ctx, i) + 2 * (hts_pos_t)ctx->half[i] + 1,
    };
}

/* Returns the cell count of row i. */
static hts_pos_t width_at(const context *ctx, size_t i)
{
    return ctx->scratch->bands[i].width;
}

/* Returns the reference prefix length a row's first cell stands for. */
static hts_pos_t origin_of(const context *ctx, size_t i)
{
    return ctx->scratch->bands[i].origin;
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

/* Returns the model's column of reference base b, which the columns cover from
 * -(margin + 1) to len + margin. */
static const phmm_column *column_of(const phmm_model *model, hts_pos_t b)
{
    return model->columns + (b + model->margin + 1);
}

/* Returns the columns of a row, such that cell k pairs column k and steps into the
 * pairing of column k + 1. prepare guarantees that every base a row's cells pair or step
 * into has a column. */
static const phmm_column *columns_of(const context *ctx, size_t i)
{
    return column_of(ctx->model, origin_of(ctx, i) - 1);
}

/* ------------------------------------------------------------------------ */
/* Forward                                                                   */
/* ------------------------------------------------------------------------ */

/* Fills the first row of the forward pass and returns its total. Cell k steps
 * only to the pairing of reference base position_of(k), so it takes that base's
 * value in ends. The alignment never begins in a deletion. The total rescales
 * the row, so the values need not sum to one. */
static double forward_first_row(const context *ctx)
{
    band_cell *row   = row_of(ctx, 0);
    hts_pos_t  width = width_at(ctx, 0);
    double     total = 0.0;

    for (hts_pos_t k = 0; k < width; k++) {
        double share = rate_at(ctx->rates->ends, ctx->ref->len,
                               position_of(ctx, 0, k));

        row[k][STATE_MATCH]     = share;
        row[k][STATE_INSERTION] = 0.0;
        row[k][STATE_DELETION]  = 0.0;

        total += share;
    }

    return total;
}

/* The row above one forward row and the transition weights that hold along the row: out
 * of the row above, with the row's factor and the insertion emission folded in, and
 * along the row itself. The decision weights vary per cell and are taken at each cell
 * instead. The folded factor is why this is the one place a forward row is read without
 * forward_at.
 *
 * The weights are copied out of the rates so that the loop over the row reads them from
 * registers; read through the rates, each would be loaded again at every cell, since
 * the compiler must assume the row's stores can reach them. The factor arrives from the
 * caller rather than from scratch for the same reason. */
typedef struct {
    const band_cell *above;
    hts_pos_t        width;
    hts_pos_t        shift;
    double           scale;
    double           match_to_insertion;
    double           insertion_to_insertion;
    double           match_to_deletion;
    double           deletion_to_deletion;
} descent;

static descent descent_into(const context *ctx, size_t i, double scale)
{
    const phmm_rates *rates = ctx->rates;

    return (descent){
        .above                  = read_row_of(ctx, i - 1),
        .width                  = width_at(ctx, i - 1),
        .shift                  = shift_between(ctx, i - 1, i),
        .scale                  = scale,
        .match_to_insertion     = rates->match_to_insertion * scale
                                * UNINFORMATIVE,
        .insertion_to_insertion = rates->insertion_to_insertion * scale
                                * UNINFORMATIVE,
        .match_to_deletion      = rates->match_to_deletion,
        .deletion_to_deletion   = rates->deletion_to_deletion,
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
static double deleted_from(const descent *step, double left_match,
                           double left_deletion)
{
    return step->match_to_deletion    * left_match
         + step->deletion_to_deletion * left_deletion;
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

/* Fills row i of the forward pass, stored with the given factor, and returns its
 * total. Valid only for a row with a row above it. */
static double forward_row(const context *ctx, size_t i, double scale)
{
    descent           step       = descent_into(ctx, i, scale);
    band_cell        *row        = row_of(ctx, i);
    cell_terms       *terms      = terms_of(ctx, i);
    const phmm_column *cols      = columns_of(ctx, i);
    row_terms         each       = row_terms_of(ctx, i);
    hts_pos_t         width      = width_at(ctx, i);
    bool              deletions  = deletions_live(ctx, i);
    bool              insertions = insertions_live(ctx, i);
    /* The cell to the left, held in locals so each step of the deletion chain
     * does not wait on the preceding store. */
    double            left_match    = 0.0;
    double            left_deletion = 0.0;
    /* One running sum per state keeps the addition chains short. */
    double            total_paired   = 0.0;
    double            total_inserted = 0.0;
    double            total_deleted  = 0.0;

    for (hts_pos_t k = 0; k < width; k++) {
        const phmm_column *col = &cols[k];  /* the base the cell pairs */
        double             paired, inserted, deleted;

        terms[k] = terms_at(&each, col);

        paired   = paired_from(&step, &col->decision, k, terms[k].emission);
        inserted = insertions ? inserted_from(&step, k) : 0.0;
        deleted  = deletions && k > 0
                 ? deleted_from(&step, left_match, left_deletion)
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

/* Checks a row's total and leaves its reciprocal in scale. A zero or subnormal total
 * means the band admits no path for this read, which is not a fault. A negative or
 * non-finite total is a fault. */
static phmm_status check_total(double total, double *scale)
{
    if (total < 0.0 || !isfinite(total)) {
        return PHMM_UNSOUND;
    }

    if (!isnormal(total)) {
        return PHMM_NO_PATH;
    }

    *scale = 1.0 / total;

    return PHMM_OK;
}

/* The reciprocals of the totals of the two rows above the one being filled, and the
 * factor of the row above, from which the factor of the next row follows. Carried in
 * registers across the rows. */
typedef struct {
    double older;   /* of the row two above */
    double newer;   /* of the row above */
    double factor;  /* the row above was stored multiplied by this */
} recent_rows;

/* Returns the factor of row i, given what is known of the rows above it. */
static double next_factor(const recent_rows *recent, size_t i)
{
    return i < 2 ? 1.0 : recent->older / recent->factor;
}

/* Fills row i of the forward pass with its factor, records the factor for the backward
 * pass, and checks the row's total. */
static phmm_status forward_step(const context *ctx, size_t i, recent_rows *recent)
{
    double      factor = next_factor(recent, i);
    double      total  = i == 0 ? forward_first_row(ctx) : forward_row(ctx, i, factor);
    double      scale  = 0.0;  /* left where the total fails, which ends the pass */
    phmm_status status = check_total(total, &scale);

    ctx->scratch->factor[i] = factor;

    recent->older  = recent->newer;
    recent->newer  = scale;
    recent->factor = factor;

    return status;
}

/* Runs the forward pass, leaving the reciprocal of the last row's total in the
 * scratch for the accumulation to restore the rows with. */
static phmm_status forward(const context *ctx)
{
    recent_rows recent = { .older = 1.0, .newer = 1.0, .factor = 1.0 };
    phmm_status status = PHMM_OK;

    for (size_t i = 0; status == PHMM_OK && i < ctx->rows; i++) {
        status = forward_step(ctx, i, &recent);
    }

    ctx->scratch->restore = recent.newer;

    return status;
}

/* ------------------------------------------------------------------------ */
/* Accumulating a row into the window                                        */
/* ------------------------------------------------------------------------ */

/* Returns the window advanced to where the row's first cell enters it, so a cell
 * addresses its positions by its own index with no bounds check. lay_out_rows
 * guarantees that every position a row can address lies inside the window. The band is
 * not clamped, so a row near either end of the reference addresses positions outside
 * it; those take only zero contributions. */
static phmm_position *landing_of(const context *ctx, size_t i)
{
    return ctx->scratch->window + (size_t)(origin_of(ctx, i) - 1 - ctx->window.origin);
}

/* The state of one row's accumulation into the window. Cells are accumulated
 * right to left, the order the backward pass forms them in. The pairing of a
 * position and the deletion of it are held by the cell to its right, so a
 * position is complete only once the cell to its left is reached. */
typedef struct {
    scaled_row        front;
    const cell_terms *terms;
    phmm_position    *at;
    /* Whether the row's pairings are also the read's 5'-most, which no decision
     * follows. */
    bool              ends;
    /* The pending contribution to the position the next cell completes. */
    double            coverage;
    double            mismatches;
    double            deletions;
} accumulation;

static accumulation accumulation_of(const context *ctx, size_t i)
{
    return (accumulation){
        .front = scaled_row_of(ctx, i),
        .terms = terms_of(ctx, i),
        .at    = landing_of(ctx, i),
        /* Row 1 pairs the read's 5'-most placed base. */
        .ends  = i == 1,
    };
}

/* Accumulates cell k of the row into the window. The two runs a cell can open begin
 * on the decision at the pairing below, whose weights the caller passes in. */
static void accumulate_cell(accumulation *acc, hts_pos_t k, const decision *below,
                            const double *back, double pairing)
{
    phmm_position *next     = &acc->at[k + 1];  /* the position this cell completes */
    double         matched  = forward_at(&acc->front, k, STATE_MATCH);
    double         skipped  = forward_at(&acc->front, k, STATE_DELETION);
    double         carried  = forward_at(&acc->front, k, STATE_INSERTION);
    double         paired   = matched * back[STATE_MATCH];
    double         deleted  = below->deletion_to_match * skipped * pairing;
    double         inserted = below->insertion_to_match * carried * pairing;

    /* Accumulated from the previous call, except the insertion. */
    next->coverage   += acc->coverage;
    next->mismatches += acc->mismatches;
    next->deletions  += acc->deletions;
    next->insertions += inserted;

    /* Read by the next call. */
    acc->coverage   = paired;
    acc->mismatches = paired * acc->terms[k].modification;
    acc->deletions  = deleted;

    /* Written at the index the pending values of this cell land at. */
    if (acc->ends) {
        acc->at[k].ends += paired;
    }
}

/* Writes the pending contribution after the leftmost cell, which completes
 * the row's last open position. */
static void accumulate_end(const accumulation *acc)
{
    acc->at[0].coverage   += acc->coverage;
    acc->at[0].mismatches += acc->mismatches;
    acc->at[0].deletions  += acc->deletions;
}

/* ------------------------------------------------------------------------ */
/* Backward                                                                  */
/* ------------------------------------------------------------------------ */

/* A backward cell holds the chance of every way the alignment can finish from that
 * cell. Each step into a row multiplies in the factor that row was stored with, so a
 * forward cell read through forward_at times a backward cell is a posterior. Each row
 * is accumulated into the window as it is formed. */

/* The row below one backward row: its cells, comparisons, the factor it was stored
 * with, its width, and the shift between the rows, with the transition weights that
 * hold along the row copied out of the rates for the reason descent gives. */
typedef struct {
    const band_cell  *cell;
    const cell_terms *terms;
    double            scale;
    hts_pos_t         width;
    hts_pos_t         shift;
    double            match_to_insertion;
    double            match_to_deletion;
    double            insertion_to_insertion;
    double            deletion_to_deletion;
} ascent;

static ascent ascent_into(const context *ctx, size_t i)
{
    const phmm_rates *rates = ctx->rates;

    return (ascent){
        .cell                   = read_backward_row_of(ctx, i + 1),
        .terms                  = terms_of(ctx, i + 1),
        .scale                  = factor_into(ctx, i + 1),
        .width                  = width_at(ctx, i + 1),
        .shift                  = shift_between(ctx, i, i + 1),
        .match_to_insertion     = rates->match_to_insertion,
        .match_to_deletion      = rates->match_to_deletion,
        .insertion_to_insertion = rates->insertion_to_insertion,
        .deletion_to_deletion   = rates->deletion_to_deletion,
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
static void backward_cell(const ascent *step, const decision *below,
                          bool insertions, double pairing, double inserted,
                          double deleted, double *cell)
{
    cell[STATE_MATCH] = below->match_to_match    * pairing
                      + step->match_to_insertion * inserted
                      + step->match_to_deletion  * deleted;

    cell[STATE_INSERTION] = insertions
                          ? below->insertion_to_match    * pairing
                          + step->insertion_to_insertion * inserted
                          : 0.0;

    cell[STATE_DELETION] = below->deletion_to_match   * pairing
                         + step->deletion_to_deletion * deleted;
}

/* Fills the last row of the backward pass and accumulates it. The alignment ends on
 * this row, and it ends on a pairing, so a match finishes with chance one and neither
 * run can be open. */
static void backward_last_row(const context *ctx)
{
    size_t             i    = ctx->rows - 1;
    band_cell         *row  = backward_row_of(ctx, i);
    const phmm_column *cols = columns_of(ctx, i);
    accumulation       acc  = accumulation_of(ctx, i);
    double             cell[N_STATES] = {
        [STATE_MATCH]     = 1.0,
        [STATE_INSERTION] = 0.0,
        [STATE_DELETION]  = 0.0,
    };

    for (hts_pos_t k = width_at(ctx, i); k-- > 0; ) {
        const decision *below = &cols[k + 1].decision;

        row[k][STATE_MATCH]     = cell[STATE_MATCH];
        row[k][STATE_INSERTION] = cell[STATE_INSERTION];
        row[k][STATE_DELETION]  = cell[STATE_DELETION];

        accumulate_cell(&acc, k, below, cell, 0.0);
    }

    accumulate_end(&acc);
}

/* Fills row i of the backward pass and accumulates it. Valid only for a row
 * with a row above and a row below. */
static void backward_row(const context *ctx, size_t i)
{
    band_cell         *row        = backward_row_of(ctx, i);
    ascent             below      = ascent_into(ctx, i);
    const phmm_column *cols       = columns_of(ctx, i);
    accumulation       acc        = accumulation_of(ctx, i);
    hts_pos_t          width      = width_at(ctx, i);
    bool               insertions = insertions_live(ctx, i);
    /* The cell to the right, held in a local so each step of the deletion
     * chain does not wait on the preceding store. */
    double            right_deletion = 0.0;

    for (hts_pos_t k = width; k-- > 0; ) {
        const decision *into     = &cols[k + 1].decision;
        double          pairing  = pairing_below(&below, k);
        double          inserted = inserted_below(&below, k);
        double          deleted  = k + 1 < width ? right_deletion : 0.0;
        double          cell[N_STATES];

        backward_cell(&below, into, insertions, pairing, inserted, deleted, cell);

        row[k][STATE_MATCH]     = cell[STATE_MATCH];
        row[k][STATE_INSERTION] = cell[STATE_INSERTION];
        row[k][STATE_DELETION]  = cell[STATE_DELETION];

        accumulate_cell(&acc, k, into, cell, pairing);

        right_deletion = cell[STATE_DELETION];
    }

    accumulate_end(&acc);
}

/* Fills the first row of the backward pass. The row precedes the first placed
 * base, so it adds no value to the window and is written only for
 * passes_agree to check. */
static void backward_first_row(const context *ctx)
{
    band_cell         *row   = backward_row_of(ctx, 0);
    ascent             below = ascent_into(ctx, 0);
    const phmm_column *cols  = columns_of(ctx, 0);
    hts_pos_t          width = width_at(ctx, 0);

    for (hts_pos_t k = 0; k < width; k++) {
        const decision *into     = &cols[k + 1].decision;
        double          pairing  = pairing_below(&below, k);
        double          inserted = inserted_below(&below, k);

        row[k][STATE_MATCH] = into->match_to_match     * pairing
                            + below.match_to_insertion * inserted;

        /* Neither run reaches the row before the first placed base. */
        row[k][STATE_INSERTION] = 0.0;
        row[k][STATE_DELETION]  = 0.0;
    }
}

/* Returns the first row's total of forward times backward, which is one when
 * the two passes describe the same set of paths. The first row carries the
 * whole backward chain and every factor, so this checks both passes end to end. */
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
/* The model                                                                 */
/* ------------------------------------------------------------------------ */

/* Returns the column of reference base b of a reference of len bases. */
static phmm_column column_at(const phmm_rates *rates, const cm_fasta_record *ref,
                             hts_pos_t b)
{
    size_t    len     = ref->len;
    bool      present = b >= 0 && (size_t)b < len;
    hts_pos_t nearest = clamp_base(len, b);

    return (phmm_column){
        .decision     = decision_at(rates, len, b),
        .modification = rates->modification[nearest],
        .code         = present ? (uint8_t)nuc_from_char(ref->seq[nearest])
                                : COLUMN_ABSENT,
    };
}

/* Returns the columns a reference of len bases needs with the given margin. */
static size_t columns_needed(size_t len, int margin)
{
    return len + 2 * (size_t)margin + 2;
}

int phmm_model_prepare(phmm_model *model, const phmm_rates *rates,
                       const cm_fasta_record *ref, int margin)
{
    size_t needed = columns_needed(ref->len, margin);

    if (needed > model->cap) {
        phmm_column *columns = realloc(model->columns, needed * sizeof *columns);

        if (!columns) {
            return -1;
        }

        model->columns = columns;
        model->cap     = needed;
    }

    model->rates  = rates;
    model->len    = ref->len;
    model->margin = margin;

    for (size_t t = 0; t < needed; t++) {
        model->columns[t] = column_at(rates, ref, (hts_pos_t)t - margin - 1);
    }

    return 0;
}

void phmm_model_free(phmm_model *model)
{
    free(model->columns);
    model->columns = NULL;
    model->cap     = 0;
}

/* ------------------------------------------------------------------------ */
/* Scratch                                                                   */
/* ------------------------------------------------------------------------ */

/* Marks every quality slot as holding no rate. */
static void clear_qualities(phmm_scratch *scratch)
{
    for (size_t q = 0; q < QUALITY_SLOTS; q++) {
        scratch->by_quality[q].modification = (double)NAN;
    }
}

phmm_scratch *phmm_scratch_create(void)
{
    phmm_scratch *scratch = calloc(1, sizeof *scratch);

    if (scratch) {
        clear_qualities(scratch);
    }

    return scratch;
}

void phmm_scratch_destroy(phmm_scratch *scratch)
{
    if (!scratch) {
        return;
    }

    free(scratch->places);
    free(scratch->bands);
    free(scratch->forward);
    free(scratch->terms);
    free(scratch->factor);
    free(scratch->backward);
    free(scratch->window);
    free(scratch);
}

/* Grows the per-row buffers. A buffer that grew is kept even when a later one
 * fails, so a failed grow leaves the scratch usable; the buffers never shrink. This
 * runs before grow_band, since the row count comes from the CIGAR while the
 * widest row is known only after the places are written. */
static int grow_rows(phmm_scratch *scratch, size_t rows)
{
    aln_place *places;
    row_band  *bands;
    double    *factor;

    if (rows <= scratch->rows) {
        return 0;
    }

    places = realloc(scratch->places, rows * sizeof *places);
    bands  = realloc(scratch->bands, rows * sizeof *bands);
    factor = realloc(scratch->factor, rows * sizeof *factor);

    if (places) {
        scratch->places = places;
    }
    if (bands) {
        scratch->bands = bands;
    }
    if (factor) {
        scratch->factor = factor;
    }

    if (!places || !bands || !factor) {
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

/* Grows the window, which holds one record per position. */
static int grow_window(phmm_scratch *scratch, size_t len)
{
    phmm_position *window;

    if (len <= scratch->window_len) {
        return 0;
    }

    window = realloc(scratch->window, len * sizeof *window);

    if (!window) {
        return -1;
    }

    scratch->window     = window;
    scratch->window_len = len;

    return 0;
}

/* ------------------------------------------------------------------------ */
/* One read                                                                  */
/* ------------------------------------------------------------------------ */

/* Writes where every row lies, and from the rows the stride every row is stored at,
 * which is the width of the widest row and at least one, and the smallest reference
 * range that holds every position any row can write. */
static void lay_out_rows(context *ctx)
{
    row_band *bands  = ctx->scratch->bands;
    hts_pos_t widest = 1;
    hts_pos_t first  = 0;
    hts_pos_t last   = 0;

    for (size_t i = 0; i < ctx->rows; i++) {
        row_band  band = band_of(ctx, i);
        hts_pos_t end  = band.origin + band.width;

        bands[i] = band;

        if (band.width > widest) {
            widest = band.width;
        }
        if (i == 0 || band.origin < first) {
            first = band.origin;
        }
        if (i == 0 || end > last) {
            last = end;
        }
    }

    ctx->widest = widest;
    ctx->window = (extent){
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

    memset(scratch->window, 0, len * sizeof *scratch->window);
}

/* Returns whether every base the read's rows reach has a column in the model. */
static bool covered(const context *ctx)
{
    const phmm_model *model = ctx->model;
    hts_pos_t         first = ctx->window.origin;
    hts_pos_t         last  = first + (hts_pos_t)ctx->window.len;

    return first >= -(model->margin + 1)
        && last <= (hts_pos_t)columns_needed(model->len, model->margin) - model->margin - 1;
}

static phmm_status prepare(context *ctx)
{
    const cm_bam_record *read = ctx->read;

    if (grow_rows(ctx->scratch, (size_t)read->l_qseq + 1) < 0) {
        return PHMM_NO_MEMORY;
    }

    ctx->span = aln_places(read, ctx->scratch->places);
    ctx->rows = (size_t)(ctx->span.end - ctx->span.begin) + 1;

    lay_out_rows(ctx);

    if (!covered(ctx)) {
        return PHMM_NO_PATH;
    }

    if (grow_band(ctx->scratch, ctx->rows, (size_t)ctx->widest) < 0) {
        return PHMM_NO_MEMORY;
    }

    if (grow_window(ctx->scratch, ctx->window.len) < 0) {
        return PHMM_NO_MEMORY;
    }

    clear_window(ctx);

    return PHMM_OK;
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

phmm_status phmm_run(const phmm_model *model, const phred *quality,
                     const cm_bam_record *read, const cm_fasta_record *ref,
                     const int *half, phmm_scratch *scratch, phmm_window *out)
{
    context ctx = {
        .model   = model,
        .rates   = model->rates,
        .quality = quality,
        .read    = read,
        .ref     = ref,
        .scratch = scratch,
        .half    = half,
    };
    phmm_status status;

    status = prepare(&ctx);

    if (status != PHMM_OK) {
        return status;
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
    out->at     = scratch->window;

    return PHMM_OK;
}
