/* cli.h -- command line parsing.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stdio.h>

#include "pipeline.h"

#define CMUTS_VERSION "0.1.0"

/* Everything the command line can set.
 *
 * Module configuration is held by nested structs rather than flattened, so
 * that each module keeps owning its own settings and this stays a container.
 * As filtering and processing options arrive they get their own structs here,
 * not more loose fields. */
typedef struct {
    pipeline_config pipeline;
    bool            show_help;
    bool            show_version;
    bool            dump_options;
} cli_args;

typedef enum {
    CLI_OK,     /* arguments parsed; carry on */
    CLI_DONE,   /* the request was answered in full; exit successfully */
    CLI_ERROR,  /* usage error, already reported */
} cli_status;

cli_args   cli_defaults(void);
cli_status cli_parse(int argc, char **argv, cli_args *args);
void       cli_usage(FILE *out, const char *program);

/* Emits the option and positional tables as JSON, for generating manual
 * pages, documentation and shell completions from the binary itself rather
 * than from a description of it kept somewhere else. */
void cli_dump_options(FILE *out);
