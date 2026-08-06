/* sample.h -- deterministic sampling from parameter specifications.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* A reproducible generator. Everything the tool decides comes from one of
 * these, so a run is a function of its seed alone. */
typedef struct {
    uint64_t state;
} rng;

void     rng_seed(rng *r, uint64_t seed);
uint64_t rng_next(rng *r);

/* Uniform on [low, high], and on [0, 1). */
long   rng_between(rng *r, long low, long high);
double rng_fraction(rng *r);

/* True with the given probability. */
bool rng_chance(rng *r, double probability);

/* How a parameter varies, written one of three ways:
 *
 *     40          every draw is 40
 *     20:200      uniform over the range, inclusive
 *     0,1,30,255  drawn from the listed values
 *
 * A constant is what pins a value down for a test; a range or a list is what
 * spreads it for a benchmark or a fuzz run. One grammar serves both, so no
 * parameter needs a separate switch to randomise it. */
#define SPEC_MAX_VALUES 64

typedef struct {
    long   low;                      /* range endpoints, equal when constant */
    long   high;
    long   values[SPEC_MAX_VALUES];  /* the list form */
    size_t n_values;                 /* zero unless the list form was used */
} spec;

/* Parses text into out. Returns 0, or -1 with a description in error. */
int spec_parse(spec *out, const char *text, char *error, size_t error_len);

/* A spec fixed to one value, for defaults. */
spec spec_constant(long value);

long spec_draw(const spec *s, rng *r);

/* The largest value a spec can produce, for sizing buffers up front. */
long spec_maximum(const spec *s);

/* Renders a spec back into its written form. */
void spec_format(const spec *s, char *out, size_t size);
