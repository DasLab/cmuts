/* params.c -- the pair HMM's rates, and the file they are read from.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "params.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>

/* Long enough for a key, a rate written out, and the spacing between them. A longer line
 * is refused, so a rate is never read from part of one. */
#define PARAMS_LINE_MAX 256

phmm_params phmm_defaults(void)
{
    return (phmm_params){
        .open_deletion    = 1e-3,
        .open_insertion   = 2e-4,
        .extend_deletion  = 0.35,
        .extend_insertion = 0.35,
        .modification     = 5e-3,
    };
}

/* Where each rate sits in the struct, by the name the file gives it. */
typedef struct {
    const char *name;
    size_t      offset;
} params_entry;

#define ENTRY(key, field) { key, offsetof(phmm_params, field) }

static const params_entry ENTRIES[] = {
    ENTRY("deletion-open",     open_deletion),
    ENTRY("deletion-extend",   extend_deletion),
    ENTRY("insertion-open",    open_insertion),
    ENTRY("insertion-extend",  extend_insertion),
    ENTRY("modification-rate", modification),
};

#define N_ENTRIES ((int)(sizeof ENTRIES / sizeof *ENTRIES))

static double *field_of(phmm_params *params, const params_entry *entry)
{
    return (double *)((char *)params + entry->offset);
}

static double value_of(const phmm_params *params, const params_entry *entry)
{
    return *(const double *)((const char *)params + entry->offset);
}

static const params_entry *entry_named(const char *name)
{
    for (int i = 0; i < N_ENTRIES; i++) {
        if (strcmp(ENTRIES[i].name, name) == 0) {
            return &ENTRIES[i];
        }
    }

    return NULL;
}

/* ------------------------------------------------------------------------ */
/* Reading                                                                   */
/* ------------------------------------------------------------------------ */

static int fail(char *error, size_t error_len, const char *path, int line,
                const char *what)
{
    snprintf(error, error_len, "%s:%d: %s", path, line, what);
    return -1;
}

/* Returns the text with the leading blanks passed over. */
static char *skip_blanks(char *text)
{
    while (*text == ' ' || *text == '\t' || *text == '\r' || *text == '\n') {
        text++;
    }

    return text;
}

/* Cuts the line at its first comment mark or line ending, and trims the blanks either
 * side of what is left. Returns the text, which is empty where the line carried none.
 *
 * The scan is bounded by the buffer so that the cut lands inside it whatever the file
 * holds. size counts the terminator, so the last index written is size - 1. */
static char *content_of(char *line, size_t size)
{
    size_t len = 0;

    while (len + 1 < size && line[len] != '\0' && line[len] != '#'
           && line[len] != '\r' && line[len] != '\n') {
        len++;
    }

    while (len > 0 && (line[len - 1] == ' ' || line[len - 1] == '\t')) {
        len--;
    }

    line[len] = '\0';

    return skip_blanks(line);
}

/* Reads one rate into params. The line holds a name and a value separated by blanks. */
static int read_entry(char *line, phmm_params *params, const char *path, int number,
                      char *error, size_t error_len)
{
    char  *name  = line;
    char  *value = line;
    char  *end   = NULL;
    double parsed;

    while (*value != '\0' && *value != ' ' && *value != '\t') {
        value++;
    }

    if (*value == '\0') {
        return fail(error, error_len, path, number, "no value for this rate");
    }

    *value = '\0';
    value  = skip_blanks(value + 1);

    const params_entry *entry = entry_named(name);

    if (!entry) {
        snprintf(error, error_len, "%s:%d: no rate is named \"%s\"", path, number, name);
        return -1;
    }

    errno  = 0;
    parsed = strtod(value, &end);

    if (*value == '\0' || !end || *end != '\0' || errno == ERANGE) {
        return fail(error, error_len, path, number, "the value is not a number");
    }

    *field_of(params, entry) = parsed;

    return 0;
}

/* Returns whether every rate is a probability, and whether the two opening rates leave
 * anything for a pairing. */
static int check(const phmm_params *params, const char *path, char *error,
                 size_t error_len)
{
    for (int i = 0; i < N_ENTRIES; i++) {
        double rate = value_of(params, &ENTRIES[i]);

        if (!(rate >= 0.0 && rate <= 1.0)) {
            snprintf(error, error_len, "%s: %s is not between 0 and 1",
                     path, ENTRIES[i].name);
            return -1;
        }
    }

    if (params->open_deletion + params->open_insertion > 1.0) {
        snprintf(error, error_len,
                 "%s: deletion-open and insertion-open sum above 1", path);
        return -1;
    }

    return 0;
}

int params_read(const char *path, phmm_params *params, char *error, size_t error_len)
{
    char        line[PARAMS_LINE_MAX];
    phmm_params read   = *params;
    int         number = 0;
    FILE       *file   = fopen(path, "r");

    if (!file) {
        /* Single-threaded, since this is reached before any thread is started.
         * NOLINTNEXTLINE(concurrency-mt-unsafe) */
        snprintf(error, error_len, "%s: %s", path, strerror(errno));
        return -1;
    }

    while (fgets(line, sizeof line, file)) {
        char *content;

        number++;

        if (!strchr(line, '\n') && !feof(file)) {
            fclose(file);
            return fail(error, error_len, path, number, "the line is too long");
        }

        content = content_of(line, sizeof line);

        if (*content == '\0') {
            continue;
        }

        if (read_entry(content, &read, path, number, error, error_len) < 0) {
            fclose(file);
            return -1;
        }
    }

    fclose(file);

    if (check(&read, path, error, error_len) < 0) {
        return -1;
    }

    *params = read;

    return 0;
}

/* ------------------------------------------------------------------------ */
/* Writing                                                                   */
/* ------------------------------------------------------------------------ */

static void params_write(const phmm_params *params, FILE *out)
{
    fprintf(out, "# The rates the pair HMM runs on, as --params reads them.\n");

    for (int i = 0; i < N_ENTRIES; i++) {
        fprintf(out, "%-18s %g\n", ENTRIES[i].name, value_of(params, &ENTRIES[i]));
    }
}

void params_dump_defaults(FILE *out)
{
    phmm_params defaults = phmm_defaults();

    params_write(&defaults, out);
}
