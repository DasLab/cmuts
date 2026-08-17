/* sample.h -- deterministic sampling from parameter specifications.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* A reproducible generator. Every value the tool draws comes from one of these, so a run is
 * a function of its seed only. */
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

/* The values a parameter may take, drawn from uniformly, written one of three ways:
 *
 *     40          every draw is 40
 *     20:200      uniform over the range, inclusive
 *     0,1,30,255  drawn from the listed values
 *
 * All three describe a finite set: a constant is a set of one, a range is the whole span
 * between its endpoints, and a list is its own contents. Only the uniform is expressible; there
 * is no way to weight one value above another.
 *
 * A constant pins a value down for a test; a range or a list spreads it for a benchmark or a
 * fuzz run. One grammar serves both, so no parameter needs a separate switch to randomize
 * it. */
#define DISTRIBUTION_MAX_VALUES 64

typedef struct {
    long   low;       /* range endpoints, equal when constant */
    long   high;
    long   values[DISTRIBUTION_MAX_VALUES];  /* the list form */
    size_t n_values;  /* zero unless the list form was used */
} distribution;

/* Parses text into out. Returns 0, or -1 with a description in error. */
int distribution_parse(distribution *out, const char *text,
                       char *error, size_t error_len);

long distribution_draw(const distribution *d, rng *r);

/* The largest value a distribution can produce, for sizing buffers up front. */
long distribution_maximum(const distribution *d);
