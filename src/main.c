/* main.c -- driver for the processing pipeline.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include <stdio.h>

#include "cli.h"
#include "pipeline.h"

#define ERROR_MAX 512

static void print_stats(const pipeline_stats *stats)
{
    printf("  reads         %zu read, %zu unmapped, %zu filtered, %zu processed\n",
           stats->reads_total, stats->reads_unmapped, stats->reads_filtered,
           stats->reads_processed);
    printf("  references    %zu written\n", stats->refs_completed);
    printf("  coverage      %.0f\n", stats->coverage_total);
    printf("  mutations     %.0f\n", stats->mutations_total);
    printf("  reads merged  %.0f\n", stats->reads_recorded);
}

int main(int argc, char **argv)
{
    cli_args       args;
    pipeline_stats stats = { 0 };
    char           error[ERROR_MAX];

    switch (cli_parse(argc, argv, &args)) {
        case CLI_DONE:  return 0;
        case CLI_ERROR: return 2;
        case CLI_OK:    break;
    }

    if (pipeline_run(&args.pipeline, &stats, error, sizeof error) < 0) {
        fprintf(stderr, "%s: %s\n", argv[0], error);
        return 1;
    }

    printf("%s -> %s, %zu workers\n",
           args.pipeline.bam_path, args.pipeline.output_path, args.pipeline.workers);
    print_stats(&stats);

    return 0;
}
