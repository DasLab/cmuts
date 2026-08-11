/* sample.c -- deterministic sampling from parameter specifications.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "sample.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------ */
/* Generator                                                                 */
/* ------------------------------------------------------------------------ */

/* splitmix64: a small generator with no state beyond a counter, which makes a
 * run reproducible across platforms without depending on the C library's. */
void rng_seed(rng *r, uint64_t seed)
{
    r->state = seed;
}

uint64_t rng_next(rng *r)
{
    uint64_t z = (r->state += 0x9e3779b97f4a7c15ULL);

    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;

    return z ^ (z >> 31);
}

/* The top bits are the ones with full quality, and this many of them are
 * exactly what a double's mantissa holds. The shift drops the rest and the
 * scale is two to the same power, so the two are written once and follow. */
#define FRACTION_BITS  53
#define FRACTION_SHIFT (64 - FRACTION_BITS)
#define FRACTION_SCALE (1.0 / (double)(1ULL << FRACTION_BITS))

double rng_fraction(rng *r)
{
    return (double)(rng_next(r) >> FRACTION_SHIFT) * FRACTION_SCALE;
}

long rng_between(rng *r, long low, long high)
{
    if (high <= low)
        return low;

    return low + (long)(rng_next(r) % (uint64_t)(high - low + 1));
}

bool rng_chance(rng *r, double probability)
{
    if (probability <= 0.0)
        return false;
    if (probability >= 1.0)
        return true;

    return rng_fraction(r) < probability;
}

/* ------------------------------------------------------------------------ */
/* Specifications                                                            */
/* ------------------------------------------------------------------------ */

/* Enough for a number as anyone would write one. A longer run of characters is
 * rejected, not truncated: a truncated run of digits parses as a different
 * number, which would be worse than failing outright. */
#define NUMBER_MAX 32

static int fail(char *error, size_t error_len, const char *text, const char *why)
{
    snprintf(error, error_len, "\"%s\" %s", text, why);
    return -1;
}

static int parse_long(const char *text, const char *end, long *out)
{
    char  buffer[NUMBER_MAX];
    char *stop = NULL;
    size_t n   = (size_t)(end - text);

    if (n == 0 || n >= sizeof buffer)
        return -1;

    memcpy(buffer, text, n);
    buffer[n] = '\0';

    *out = strtol(buffer, &stop, 10);
    return *stop == '\0' ? 0 : -1;
}

static int parse_list(distribution *out, const char *text,
                      char *error, size_t error_len)
{
    const char *from = text;

    out->n_values = 0;

    for (;;) {
        const char *comma = strchr(from, ',');
        const char *end   = comma ? comma : from + strlen(from);

        if (out->n_values == DISTRIBUTION_MAX_VALUES)
            return fail(error, error_len, text, "lists more values than are allowed");

        if (parse_long(from, end, &out->values[out->n_values]) < 0)
            return fail(error, error_len, text, "is not a list of numbers");

        out->n_values++;

        if (!comma)
            break;
        from = comma + 1;
    }

    return 0;
}

static int parse_range(distribution *out, const char *text, const char *colon,
                       char *error, size_t error_len)
{
    if (parse_long(text, colon, &out->low) < 0 ||
        parse_long(colon + 1, text + strlen(text), &out->high) < 0)
        return fail(error, error_len, text, "is not a range of the form LOW:HIGH");

    if (out->high < out->low)
        return fail(error, error_len, text, "has its range the wrong way round");

    return 0;
}

int distribution_parse(distribution *out, const char *text,
                       char *error, size_t error_len)
{
    const char *colon = strchr(text, ':');

    *out = (distribution){ 0 };

    if (text[0] == '\0')
        return fail(error, error_len, text, "is empty");

    if (strchr(text, ','))
        return parse_list(out, text, error, error_len);

    if (colon)
        return parse_range(out, text, colon, error, error_len);

    if (parse_long(text, text + strlen(text), &out->low) < 0)
        return fail(error, error_len, text, "is not a number");

    out->high = out->low;
    return 0;
}

long distribution_draw(const distribution *d, rng *r)
{
    if (d->n_values)
        return d->values[rng_next(r) % d->n_values];

    return rng_between(r, d->low, d->high);
}

long distribution_maximum(const distribution *d)
{
    long largest;

    if (!d->n_values)
        return d->high;

    largest = d->values[0];
    for (size_t i = 1; i < d->n_values; i++)
        if (d->values[i] > largest)
            largest = d->values[i];

    return largest;
}
