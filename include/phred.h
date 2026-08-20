/* phred.h -- base qualities as probabilities.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdint.h>

/* A score arrives as one byte, so every value a byte can hold has an entry. */
#define PHRED_MAX 255

/* A table does not change after it is built, so every thread can share one. */
typedef struct {
    double error[PHRED_MAX + 1];
} phred;

/* Builds the table. A base that scores below the minimum gives no information about the
 * template. A minimum of 0 keeps every base at its own score. */
void phred_build(phred *table, int minimum);

/* Returns the chance a base call is wrong. */
static inline double phred_error(const phred *table, uint8_t quality)
{
    return table->error[quality];
}

/* Returns the chance a base call is right. */
static inline double phred_correct(const phred *table, uint8_t quality)
{
    return 1.0 - table->error[quality];
}
