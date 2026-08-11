/* progress.h -- how far through the input the reader has got.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include "bamstream.h"

typedef struct progress progress;

/* Returns NULL where there is no terminal to draw on, or where the size of the input is
 * unknown. The other two accept NULL, so a caller never has to check which. */
progress *progress_start(const cm_bam_stream *stream);

void progress_follow(progress *bar);
void progress_finish(progress *bar);
