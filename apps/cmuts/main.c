/* main.c -- dispatch to a subcommand.
 *
 * align is a script installed alongside the binary, so it is run through exec
 * rather than called.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "subcommands.h"
#include "version.h"

#define PROGRAM "cmuts"
#define ALIGN   "align"

typedef struct {
    const char *name;
    int       (*run)(int argc, char **argv);
    const char *summary;
} subcommand;

/* In pipeline order, which is the order the help lists them in. align has no
 * entry point here and is dispatched by name. */
static const subcommand SUBCOMMANDS[] = {
    { ALIGN,  NULL,      "align reads to a reference and sort the resulting alignments" },
    { "hmm",  hmm_main,  "count MaP-seq mutations via the pair-HMM" },
    { "sub",  sub_main,  "subtract an untreated background from an output" },
    { "div",  div_main,  "divide an output by a denatured control" },
    { "norm", norm_main, "divide outputs by one scale taken from their rates" },
    { "gen",  gen_main,  "generate alignments and the reference they came from" },
};

#define N_SUBCOMMANDS (sizeof SUBCOMMANDS / sizeof *SUBCOMMANDS)

static void usage(FILE *out)
{
    fprintf(out, "%s %s -- Fast, multithreaded pair-HMM counting of MaP-seq mutations.\n\n",
            PROGRAM, CMUTS_VERSION);
    fprintf(out, "usage: %s SUBCOMMAND [ARGS]\n\nSubcommands:\n", PROGRAM);

    for (size_t i = 0; i < N_SUBCOMMANDS; i++) {
        fprintf(out, "  %-6s %s\n", SUBCOMMANDS[i].name, SUBCOMMANDS[i].summary);
    }

    fprintf(out, "\nRun '%s SUBCOMMAND --help' for that subcommand's options.\n", PROGRAM);
}

/* Replaces this process with the align script, which must be on PATH. */
static int exec_align(char **argv)
{
    static char script[] = "cmuts-align";

    argv[0] = script;
    execvp(script, argv);
    /* Single-threaded, since no subcommand has run.
     * NOLINTNEXTLINE(concurrency-mt-unsafe) */
    fprintf(stderr, "%s: %s: %s\n", PROGRAM, script, strerror(errno));

    return 127;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        usage(stderr);
        return 2;
    }

    const char *name = argv[1];

    if (strcmp(name, "-h") == 0 || strcmp(name, "--help") == 0) {
        usage(stdout);
        return 0;
    }

    if (strcmp(name, "-V") == 0 || strcmp(name, "--version") == 0) {
        printf("%s %s\n", PROGRAM, CMUTS_VERSION);
        return 0;
    }

    if (strcmp(name, ALIGN) == 0) {
        return exec_align(argv + 1);
    }

    for (size_t i = 0; i < N_SUBCOMMANDS; i++) {
        if (SUBCOMMANDS[i].run && strcmp(name, SUBCOMMANDS[i].name) == 0) {
            return SUBCOMMANDS[i].run(argc - 1, argv + 1);
        }
    }

    fprintf(stderr, "%s: no subcommand named '%s'; run '%s --help' for the list\n",
            PROGRAM, name, PROGRAM);

    return 2;
}
