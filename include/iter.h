/* iter.h -- iteration protocol shared by every record reader.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

/* Return value of every reader's next() function. A reader that has returned CM_ITER_EOF or
 * CM_ITER_ERROR must not be advanced again. */
typedef enum {
    CM_ITER_OK    =  1,
    CM_ITER_EOF   =  0,
    CM_ITER_ERROR = -1,
} cm_iter_status;
