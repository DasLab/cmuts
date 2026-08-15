/* divide.c -- the rule for each field of a normalization.
 *
 * The reading and writing are combine's. This file holds the arithmetic and the table of
 * which rule applies to which field.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "divide.h"

#include <math.h>

#include "combine.h"
#include "output.h"

/* The name written into the output as the program that produced it. */
#define DIVIDE_PROGRAM "cmuts-div"

typedef enum {
    DIV_RATES,
    DIV_CONTROL,
    DIV_N_INPUTS,
} div_input;

typedef enum {
    DIV_SUM,
    DIV_RATIO,
    DIV_RATIO_ERROR,
} div_rule;

static const div_rule RULES[OUT_N_FIELDS] = {
    [OUT_COVERAGE]   = DIV_SUM,
    [OUT_REACTIVITY] = DIV_RATIO,
    [OUT_ERROR]      = DIV_RATIO_ERROR,
    [OUT_LENGTHS]    = DIV_SUM,
    [OUT_READS]      = DIV_SUM,
    [OUT_REJECTED]   = DIV_SUM,
    [OUT_UNMAPPED]   = DIV_SUM,
};

/* ------------------------------------------------------------------------ */
/* Arithmetic                                                                */
/* ------------------------------------------------------------------------ */

/* NaN in either input carries through every rule, so an unmeasured position and the
 * columns past a reference stay NaN. The result is NaN wherever the control is at or
 * below zero. */

static void ratio_f32(const float *rates, const float *control, float *out, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        out[i] = control[i] > 0.0F ? rates[i] / control[i] : (float)NAN;
    }
}

/* Returns the error of the ratio, treating the two runs as independent.
 * Not the relative-error form: that divides by the rate, which is zero at an unreactive
 * position, and squaring the control into a denominator underflows a float long before the
 * control does. */
static void ratio_error_f32(const combine_rows *rows, float *out, size_t n)
{
    const float *rate          = combine_row(rows, DIV_RATES, OUT_REACTIVITY);
    const float *control       = combine_row(rows, DIV_CONTROL, OUT_REACTIVITY);
    const float *rate_error    = combine_row(rows, DIV_RATES, OUT_ERROR);
    const float *control_error = combine_row(rows, DIV_CONTROL, OUT_ERROR);

    for (size_t i = 0; i < n; i++) {
        float d = control[i];
        float r = rate[i] / d;
        float e = rate_error[i] * rate_error[i];
        float c = r * control_error[i];

        out[i] = d > 0.0F ? sqrtf(e + c * c) / d : (float)NAN;
    }
}

/* ------------------------------------------------------------------------ */
/* Rules                                                                     */
/* ------------------------------------------------------------------------ */

static int combine_f32(const combine_rows *rows, out_field_id id, div_rule how, float *out,
                       size_t n)
{
    switch (how) {
        case DIV_RATIO:
            ratio_f32(combine_row(rows, DIV_RATES, id),
                      combine_row(rows, DIV_CONTROL, id), out, n);
            return 0;
        case DIV_RATIO_ERROR:
            ratio_error_f32(rows, out, n);
            return 0;
        case DIV_SUM:
            break;
    }

    return -1;
}

static int divide_field(const combine_rows *rows, out_field_id id, void *out, size_t n,
                        const void *ctx)
{
    (void)ctx;

    /* Every count sums; every other rule is over rates. */
    if (RULES[id] == DIV_SUM) {
        return combine_sum(rows, id, out, n);
    }

    if (OUT_FIELDS[id].stored != OUT_F32) {
        return -1;
    }

    return combine_f32(rows, id, RULES[id], out, n);
}

/* ------------------------------------------------------------------------ */
/* Assembly                                                                  */
/* ------------------------------------------------------------------------ */

int divide_run(const divide_config *cfg, char *error, size_t error_len)
{
    const char *paths[DIV_N_INPUTS];

    paths[DIV_RATES]   = cfg->rates_path;
    paths[DIV_CONTROL] = cfg->control_path;

    const combine_spec spec = {
        .program   = DIVIDE_PROGRAM,
        .inputs    = paths,
        .n_inputs  = DIV_N_INPUTS,
        .output    = cfg->output_path,
        .overwrite = cfg->overwrite,
        .field     = divide_field,
    };

    return combine_run(&spec, error, error_len);
}
