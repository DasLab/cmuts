/* subtract.c -- the rule for each field of a background subtraction.
 *
 * The reading and writing are combine's. This file holds the arithmetic and the table of
 * which rule applies to which field.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "subtract.h"

#include <math.h>
#include <stdbool.h>

#include "combine.h"
#include "output.h"

typedef enum {
    SUB_TREATED,
    SUB_UNTREATED,
    SUB_N_INPUTS,
} sub_input;

typedef enum {
    SUB_SUM,
    SUB_SAME,
    SUB_DIFFERENCE,
    SUB_QUADRATURE,
} sub_rule;

static const sub_rule RULES[OUT_N_FIELDS] = {
    [OUT_COVERAGE]   = SUB_SUM,
    [OUT_REACTIVITY] = SUB_DIFFERENCE,
    [OUT_ERROR]      = SUB_QUADRATURE,
    [OUT_LENGTHS]    = SUB_SUM,
    [OUT_READS]      = SUB_SUM,
    [OUT_REJECTED]   = SUB_SUM,
    [OUT_UNMAPPED]   = SUB_SUM,
    [OUT_SEQUENCE]   = SUB_SAME,
};

/* ------------------------------------------------------------------------ */
/* Arithmetic                                                                */
/* ------------------------------------------------------------------------ */

/* NaN in either input carries through every rule, so an unmeasured position and the
 * columns past a reference stay NaN. */

static void subtract_f32(const float *treated, const float *untreated, float *out,
                         size_t n)
{
    for (size_t i = 0; i < n; i++) {
        out[i] = treated[i] - untreated[i];
    }
}

/* Raises every negative value to zero, leaving NaN as it is. */
static void clip_f32(float *out, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        if (out[i] < 0.0F) {
            out[i] = 0.0F;
        }
    }
}

/* Each square is rounded before the two are added, which a single expression would let the
 * compiler fuse into a multiply-add. Fusing keeps one square at a wider precision than the
 * other, and the quadrature then depends on which input came first. */
static void propagate_f32(const float *treated, const float *untreated, float *out,
                          size_t n)
{
    for (size_t i = 0; i < n; i++) {
        float t = treated[i] * treated[i];
        float u = untreated[i] * untreated[i];

        out[i] = sqrtf(t + u);
    }
}

/* ------------------------------------------------------------------------ */
/* Rules                                                                     */
/* ------------------------------------------------------------------------ */

static int combine_f32(const combine_rows *rows, out_field_id id, sub_rule how, bool clip,
                       float *out, size_t n)
{
    const float *treated   = combine_row(rows, SUB_TREATED, id);
    const float *untreated = combine_row(rows, SUB_UNTREATED, id);

    switch (how) {
        case SUB_DIFFERENCE:
            subtract_f32(treated, untreated, out, n);
            if (clip) {
                clip_f32(out, n);
            }
            return 0;
        case SUB_QUADRATURE:
            propagate_f32(treated, untreated, out, n);
            return 0;
        case SUB_SUM:
        case SUB_SAME:
            break;
    }

    return COMBINE_NO_RULE;
}

static int subtract_field(const combine_rows *rows, out_field_id id, void *out, size_t n,
                          const void *ctx)
{
    const subtract_config *cfg = ctx;

    /* Every count sums and every other rule is over rates. */
    if (RULES[id] == SUB_SUM) {
        return combine_sum(rows, id, out, n);
    }

    if (RULES[id] == SUB_SAME) {
        return combine_same(rows, id, out, n);
    }

    if (OUT_FIELDS[id].stored != OUT_F32) {
        return COMBINE_NO_RULE;
    }

    return combine_f32(rows, id, RULES[id], cfg->clip, out, n);
}

/* ------------------------------------------------------------------------ */
/* Assembly                                                                  */
/* ------------------------------------------------------------------------ */

int subtract_run(const subtract_config *cfg, const char *program,
                 const out_manifest *writes, char *error,
                 size_t error_len)
{
    const char *paths[SUB_N_INPUTS];

    paths[SUB_TREATED]   = cfg->treated_path;
    paths[SUB_UNTREATED] = cfg->untreated_path;

    const combine_spec spec = {
        .program   = program,
        .writes    = writes,
        .inputs    = paths,
        .n_inputs  = SUB_N_INPUTS,
        .output    = cfg->output_path,
        .overwrite = cfg->overwrite,
        .field     = subtract_field,
        .ctx       = cfg,
    };

    return combine_run(&spec, error, error_len);
}
