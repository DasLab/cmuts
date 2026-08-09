/* phred.h -- base qualities as probabilities.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdint.h>

/* A score arrives as one byte, so every value one can hold has an entry. */
#define PHRED_MAX 255

/* Nothing writes to it once built, so one may be shared by every thread. Held
 * by whoever needs it rather than kept as a global, which is what makes that
 * true of it by construction. */
typedef struct {
    double error[PHRED_MAX + 1];
} phred;

void phred_build(phred *table);

/* The chance a base call is wrong, and the chance it is right. Read from the
 * table rather than computed, since it is needed for every base of every
 * read. */
static inline double phred_error(const phred *table, uint8_t quality)
{
    return table->error[quality];
}

static inline double phred_correct(const phred *table, uint8_t quality)
{
    return 1.0 - table->error[quality];
}
