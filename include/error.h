/* error.h -- the size of every error message buffer.
 *
 * One size for every buffer, since a message written into one is routinely copied into
 * another.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#define CM_ERROR_MAX 512
