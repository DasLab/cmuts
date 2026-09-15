/* proc.h -- child processes and the pipes between them.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <sys/types.h>

/* Pass this instead of a file descriptor to let the child keep the descriptor that this
 * process uses. */
#define PROC_INHERIT (-1)

/* Returns true if a directory in PATH holds the program and the program is executable. */
bool proc_on_path(const char *program);

/* Opens a pipe. The read end goes to fds[0] and the write end to fds[1].
 *
 * Both ends close when a child runs its program, so a child started later cannot hold an
 * end that belongs to an earlier stage. An end that proc_spawn gives to a child survives
 * in that child. Returns 0, or -1 with a description in error. */
int proc_pipe(int fds[2], char *error, size_t error_len);

/* Starts argv[0], which is searched for in PATH, and writes its process id to pid. The
 * child reads stdin_fd as its standard input and writes stdout_fd as its standard output,
 * and it starts with the default handling of every signal.
 *
 * This function closes each descriptor that is not PROC_INHERIT, because the child holds
 * it from now on. The caller therefore keeps only the pipe ends that it uses itself.
 * Returns 0, or -1 with a description in error. */
int proc_spawn(const char *const argv[], int stdin_fd, int stdout_fd, pid_t *pid,
               char *error, size_t error_len);

/* Waits for a child and reports how it ended. Returns 0 if it exited successfully, or -1
 * with a description in error that names program. */
int proc_wait(pid_t pid, const char *program, char *error, size_t error_len);
