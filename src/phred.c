/* phred.c -- the score-to-probability table.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "phred.h"

#include <math.h>

/* A PHRED score is ten times the negative log of the chance the call is wrong,
 * so a score of 0 means no confidence at all and 30 means one in a thousand. */
void phred_build(phred *table)
{
    for (unsigned quality = 0; quality <= PHRED_MAX; quality++)
        table->error[quality] = pow(10.0, -(double)quality / 10.0);
}
