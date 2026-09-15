/* proc.c -- starting child processes, and reading how they ended.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "proc.h"

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <signal.h>
#include <spawn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

extern char **environ;

/* The functions here run before the pipeline starts its threads. The library calls that
 * are not thread safe are therefore safe to use.
 * NOLINTBEGIN(concurrency-mt-unsafe) */

static int fail(const char *what, int code, char *error, size_t error_len)
{
    snprintf(error, error_len, "%s: %s", what, strerror(code));

    return -1;
}

/* ------------------------------------------------------------------------ */
/* Finding a program                                                         */
/* ------------------------------------------------------------------------ */

/* Returns true if this directory holds the program as a regular file that may be run.
 *
 * A directory is executable when it may be searched, so the execute bit alone does not
 * tell an ordinary directory apart from a program of the same name. An empty entry in
 * PATH means the current directory, which is how the shell reads it. */
static bool executable_in(const char *directory, size_t len, const char *program)
{
    char        path[PATH_MAX];
    struct stat info;

    if (len == 0) {
        directory = ".";
        len       = 1;
    }

    if (snprintf(path, sizeof path, "%.*s/%s", (int)len, directory, program)
        >= (int)sizeof path) {
        return false;
    }

    return stat(path, &info) == 0 && S_ISREG(info.st_mode) && access(path, X_OK) == 0;
}

bool proc_on_path(const char *program)
{
    const char *path = getenv("PATH");

    if (!path) {
        return false;
    }

    while (true) {
        const char *separator = strchr(path, ':');
        size_t      len       = separator ? (size_t)(separator - path) : strlen(path);

        if (executable_in(path, len, program)) {
            return true;
        }

        if (!separator) {
            return false;
        }

        path = separator + 1;
    }
}

/* ------------------------------------------------------------------------ */
/* Starting and reaping                                                      */
/* ------------------------------------------------------------------------ */

/* Both ends of a pipe are closed when a child runs its program, so that a child started
 * later does not hold an end that belongs to an earlier stage. A pipe end that a child is
 * meant to keep survives, because duplicating a descriptor clears the flag. */
static int close_on_exec(int fd)
{
    return fcntl(fd, F_SETFD, FD_CLOEXEC);
}

int proc_pipe(int fds[2], char *error, size_t error_len)
{
    if (pipe(fds) < 0) {
        return fail("pipe", errno, error, error_len);
    }

    if (close_on_exec(fds[0]) < 0 || close_on_exec(fds[1]) < 0) {
        int code = errno;

        close(fds[0]);
        close(fds[1]);

        return fail("pipe", code, error, error_len);
    }

    return 0;
}

/* Makes the child use fd as the given standard stream. */
static int redirect(posix_spawn_file_actions_t *actions, int fd, int stream)
{
    return fd == PROC_INHERIT ? 0 : posix_spawn_file_actions_adddup2(actions, fd, stream);
}

/* Closes the descriptors that the child holds from now on. A read from a pipe returns
 * the end of the stream only after every other process has closed the write end. */
static void release(int stdin_fd, int stdout_fd)
{
    if (stdin_fd != PROC_INHERIT) {
        close(stdin_fd);
    }

    if (stdout_fd != PROC_INHERIT) {
        close(stdout_fd);
    }
}

/* posix_spawnp declares its argument vector without const, but it does not write to it.
 * The union passes the same pointers without a cast that discards const. */
static char *const *unqualified(const char *const argv[])
{
    union {
        const char *const *given;
        char *const       *taken;
    } vector = { .given = argv };

    return vector.taken;
}

/* Gives the child the default handling of every signal. The caller ignores some of them
 * for its own sake, and a child that inherited that would report a write error where it
 * should simply end. */
static int default_signals(posix_spawnattr_t *attributes)
{
    sigset_t every;
    int      code = posix_spawnattr_init(attributes);

    if (code != 0) {
        return code;
    }

    sigfillset(&every);

    code = posix_spawnattr_setsigdefault(attributes, &every);

    if (code == 0) {
        code = posix_spawnattr_setflags(attributes, POSIX_SPAWN_SETSIGDEF);
    }

    return code;
}

int proc_spawn(const char *const argv[], int stdin_fd, int stdout_fd, pid_t *pid,
               char *error, size_t error_len)
{
    posix_spawn_file_actions_t actions;
    posix_spawnattr_t          attributes;
    int                        code = posix_spawn_file_actions_init(&actions);

    if (code != 0) {
        return fail(argv[0], code, error, error_len);
    }

    code = default_signals(&attributes);

    if (code == 0) {
        code = redirect(&actions, stdin_fd, STDIN_FILENO);
    }

    if (code == 0) {
        code = redirect(&actions, stdout_fd, STDOUT_FILENO);
    }

    if (code == 0) {
        code = posix_spawnp(pid, argv[0], &actions, &attributes, unqualified(argv),
                            environ);
    }

    posix_spawnattr_destroy(&attributes);
    posix_spawn_file_actions_destroy(&actions);
    release(stdin_fd, stdout_fd);

    return code == 0 ? 0 : fail(argv[0], code, error, error_len);
}

int proc_wait(pid_t pid, const char *program, char *error, size_t error_len)
{
    int status = 0;

    while (waitpid(pid, &status, 0) < 0) {
        if (errno != EINTR) {
            return fail(program, errno, error, error_len);
        }
    }

    if (WIFSIGNALED(status)) {
        snprintf(error, error_len, "%s: killed by signal %d", program, WTERMSIG(status));
        return -1;
    }

    if (WEXITSTATUS(status) != 0) {
        snprintf(error, error_len, "%s: exited with status %d", program,
                 WEXITSTATUS(status));
        return -1;
    }

    return 0;
}

/* NOLINTEND(concurrency-mt-unsafe) */
