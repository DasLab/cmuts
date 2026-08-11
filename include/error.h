/* error.h -- how much room a failure gets to describe itself.
 *
 * One size for every such buffer, a message written into one being routinely copied into
 * another: a reason from the reference reader becomes part of the pipeline's message, which
 * becomes what the program prints. Sizes that drifted apart would truncate at the seam.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#define CM_ERROR_MAX 512
