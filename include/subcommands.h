/* subcommands.h -- the entry point of each subcommand.
 *
 * Each takes the arguments from its own name onward, so argv[0] is the
 * subcommand and argv[1] its first argument.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

int hmm_main(int argc, char **argv);
int gen_main(int argc, char **argv);
int sub_main(int argc, char **argv);
int div_main(int argc, char **argv);
int norm_main(int argc, char **argv);
int score_main(int argc, char **argv);
