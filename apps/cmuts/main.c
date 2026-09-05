/* main.c -- dispatch to a subcommand.
 *
 * align and plot are scripts installed alongside the binary, so they are run
 * through exec rather than called.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "format.h"
#include "subcommands.h"
#include "version.h"

#define PROGRAM "cmuts"

/* The longest name a script subcommand resolves to: "cmuts-" and the name. */
#define SCRIPT_NAME_MAX 64

typedef struct {
    const char *name;
    int       (*run)(int argc, char **argv);
    const char *summary;
} subcommand;

/* In pipeline order, which is the order the help lists them in. An entry with
 * no entry point is a script installed alongside the binary, dispatched by
 * name as cmuts-NAME. */
static const subcommand SUBCOMMANDS[] = {
    { "align", NULL,       "align reads to a reference and sort the resulting alignments" },
    { "hmm",   hmm_main,   "count MaP-seq mutations via the pair HMM" },
    { "sub",   sub_main,   "subtract an untreated background from an output" },
    { "div",   div_main,   "divide an output by a denatured control" },
    { "norm",  norm_main,  "normalize reactivity values across experiments" },
    { "score", score_main, "measure an output against a known structure" },
    { "plot",  NULL,       "serve an interactive report over outputs" },
    { "gen",   gen_main,   "generate alignments and the reference they came from" },
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

/* Replaces this process with the named script, which must be on PATH. */
static int exec_script(const char *name, char **argv)
{
    static char script[SCRIPT_NAME_MAX];

    snprintf(script, sizeof script, "%s-%s", PROGRAM, name);
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

    /* Hidden, like the subcommands' dump flags: the format belongs to the suite rather
     * than to one program, so the dispatcher describes it. */
    if (strcmp(name, "--dump-layout") == 0) {
        fmt_dump_format(stdout);
        return 0;
    }

    for (size_t i = 0; i < N_SUBCOMMANDS; i++) {
        if (strcmp(name, SUBCOMMANDS[i].name) != 0) {
            continue;
        }

        if (SUBCOMMANDS[i].run) {
            return SUBCOMMANDS[i].run(argc - 1, argv + 1);
        }

        return exec_script(SUBCOMMANDS[i].name, argv + 1);
    }

    fprintf(stderr, "%s: no subcommand named '%s'; run '%s --help' for the list\n",
            PROGRAM, name, PROGRAM);

    return 2;
}
