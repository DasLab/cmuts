/* progress.c -- a bar on the terminal, and nothing anywhere else.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "progress.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <unistd.h>

/* Columns spent on "[", "]", a space and "100%". */
#define BAR_FURNITURE 7

#define BAR_MIN_CELLS   10
#define BAR_MAX_CELLS   60
#define ASSUMED_COLUMNS 80

/* Cells are UTF-8, so a cell occupies one column but more than one byte. */
#define CELL_FILLED "█"
#define CELL_EMPTY  " "
#define CELL_MAX_BYTES 4

struct progress {
    const cm_bam_stream *stream;
    uint64_t             span;
    uint64_t             next;   /* position at which the bar would change */
    int                  width;  /* cells, not columns */
    int                  shown;  /* percentage drawn, or -1 */
};

static int terminal_cells(void)
{
    struct winsize size;
    int            columns = ASSUMED_COLUMNS;
    int            cells;

    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &size) == 0 && size.ws_col > BAR_FURNITURE)
        columns = size.ws_col;

    cells = columns - BAR_FURNITURE;

    if (cells < BAR_MIN_CELLS)
        return BAR_MIN_CELLS;
    if (cells > BAR_MAX_CELLS)
        return BAR_MAX_CELLS;

    return cells;
}

progress *progress_start(const cm_bam_stream *stream)
{
    uint64_t  span = cm_bam_stream_span(stream);
    progress *bar;

    if (span == 0 || !isatty(STDOUT_FILENO))
        return NULL;

    bar = calloc(1, sizeof *bar);
    if (!bar)
        return NULL;

    bar->stream = stream;
    bar->span   = span;
    bar->width  = terminal_cells();
    bar->shown  = -1;

    return bar;
}

/* Decoding runs ahead of what has been handed on, so the reader can be past the
 * end of the span the bar was drawn across. */
static int percentage(const progress *bar, uint64_t position)
{
    return position >= bar->span ? 100 : (int)(position * 100 / bar->span);
}

static uint64_t redraw_at(const progress *bar, int percent)
{
    return percent >= 100 ? UINT64_MAX : bar->span * (uint64_t)(percent + 1) / 100;
}

static void draw(const progress *bar, int percent)
{
    char   cells[BAR_MAX_CELLS * CELL_MAX_BYTES + 1];
    int    filled = bar->width * percent / 100;
    size_t at     = 0;

    for (int cell = 0; cell < bar->width; cell++) {
        const char *glyph = cell < filled ? CELL_FILLED : CELL_EMPTY;
        size_t      bytes = strlen(glyph);

        memcpy(cells + at, glyph, bytes);
        at += bytes;
    }

    cells[at] = '\0';

    printf("\r[%s] %3d%%", cells, percent);
    fflush(stdout);
}

/* Asking costs a comparison until the moment there is something new to show,
 * which is why the loader may ask on every read. */
void progress_follow(progress *bar)
{
    uint64_t position;
    int      percent;

    if (!bar)
        return;

    position = cm_bam_stream_position(bar->stream);

    if (position < bar->next)
        return;

    percent = percentage(bar, position);

    if (percent != bar->shown) {
        draw(bar, percent);
        bar->shown = percent;
    }

    bar->next = redraw_at(bar, percent);
}

/* The bar stays as it ended, showing how far a failed run got. Whatever is
 * reported next begins on a line of its own. */
void progress_finish(progress *bar)
{
    if (!bar)
        return;

    if (bar->shown >= 0) {
        putchar('\n');
        fflush(stdout);
    }

    free(bar);
}
