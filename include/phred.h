/* phred.h -- base qualities as probabilities.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdint.h>

/* A score arrives as one byte, so every value a byte can hold has an entry. */
#define PHRED_MAX 255

/* A table is read-only once built, so one may be shared by every thread. Held by whoever
 * needs it rather than kept as a global. */
typedef struct {
    double error[PHRED_MAX + 1];
} phred;

void phred_build(phred *table);

/* The chance a base call is wrong, and the chance it is right. Read from the table rather than
 * computed, being needed for every base of every read. */
static inline double phred_error(const phred *table, uint8_t quality)
{
    return table->error[quality];
}

static inline double phred_correct(const phred *table, uint8_t quality)
{
    return 1.0 - table->error[quality];
}
