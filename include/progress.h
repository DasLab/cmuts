/* progress.h -- how far through its input a run has got.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdint.h>

typedef struct progress progress;

/* Starts a bar over span units of work, in whatever unit the caller measures its
 * positions. Returns NULL where there is no terminal to draw on, or where span is 0.
 * The other two accept NULL, so a caller never has to check which. */
progress *progress_start(uint64_t span);

void progress_follow(progress *bar, uint64_t position);
void progress_finish(progress *bar);
