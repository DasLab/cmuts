/* combine.h -- several outputs read together and written as one.
 *
 * Holds everything but the arithmetic: opening the inputs, checking their shapes match,
 * reading a row at a time, and writing the result. The caller supplies the rule for each
 * field.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>

#include "output.h"

/* One reference's values, one row per input per field. A field with no row holds one value
 * here, so the same rule applies to a field with a row and to one without. */
typedef struct combine_rows combine_rows;

/* Gives one input's values for a field, in the type the field is stored as. Any field of
 * any input is available, not only the one being formed. */
const void *combine_row(const combine_rows *rows, size_t input, out_field_id id);

/* Sums one field across every input. Returns 0, or -1 for a type with no sum. */
int combine_sum(const combine_rows *rows, out_field_id id, void *out, size_t n);

/* Fills out with the n values one field of the result holds. Returns 0, or -1 where no
 * rule applies to the field. */
typedef int (*combine_fn)(const combine_rows *rows, out_field_id id, void *out, size_t n,
                          const void *ctx);

typedef struct {
    const char        *program;  /* the name written into the output */
    const char *const *inputs;   /* paths, in the order the rules index them */
    size_t             n_inputs;
    const char        *output;
    bool               overwrite;
    combine_fn         field;
    const void        *ctx;      /* passed to field; holds the caller's settings */
} combine_spec;

/* Combines the inputs field by field into a file of the same layout, under the rule given
 * for each field. Returns 0, or -1 with a description in error. */
int combine_run(const combine_spec *spec, char *error, size_t error_len);
