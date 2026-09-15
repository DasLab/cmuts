/* mapping.c -- the alignment run, from the reads to a sorted output file.
 *
 * The reads reach minimap2 directly, or through fastp when they arrive as a pair of
 * files, or through this process when they arrive as a BAM file. The alignments that
 * minimap2 writes pass through this process on their way to samtools sort, which is where
 * the reference checksums and the read names are corrected.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "mapping.h"

#include <dirent.h>
#include <errno.h>
#include <libgen.h>
#include <limits.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include <htslib/kstring.h>
#include <htslib/sam.h>

#include "checksum.h"
#include "error.h"
#include "fasta.h"
#include "iter.h"
#include "proc.h"
#include "reads.h"

/* A BAM file stores the length of a read name in one byte, so the name and the NUL that
 * ends it must together fit in 255 bytes. */
#define NAME_LIMIT 254

/* fastp accepts no more threads than this. */
#define MERGE_THREADS_MAX 16

/* Room for the decimal form of a thread count. */
#define THREADS_TEXT_MAX 16

/* Room for the decimal form of a memory limit and the suffix that follows it. */
#define MEMORY_TEXT_MAX 24

/* The buffer that every line of the alignments is read into. It grows where a line does
 * not fit, so this is a starting point and not a limit. */
#define LINE_SIZE 4096

/* The first bytes of a minimap2 index. */
static const char INDEX_MAGIC[] = { 'M', 'M', 'I', '\2' };

/* How many files of reads the run was given. The command line leaves the second path
 * NULL where only one file was given. */
static size_t reads_given(const mapping_config *cfg)
{
    return cfg->reads_paths[1] ? MAPPING_MAX_READS : 1;
}

static int fail(const char *what, const char *why, char *error, size_t error_len)
{
    snprintf(error, error_len, "%s: %s", what, why);

    return -1;
}

/* Reports a failure that names no file. */
static int fail_plainly(const char *why, char *error, size_t error_len)
{
    snprintf(error, error_len, "%s", why);

    return -1;
}

/* Closes a descriptor. Where a stream was opened over it, closes the stream instead. */
static void close_either(FILE *stream, int fd)
{
    if (stream) {
        fclose(stream);
    } else {
        close(fd);
    }
}

static int fail_with_errno(const char *what, char *error, size_t error_len)
{
    /* The checks run before the pipeline starts its threads.
     * NOLINTNEXTLINE(concurrency-mt-unsafe) */
    return fail(what, strerror(errno), error, error_len);
}

/* ------------------------------------------------------------------------ */
/* Checking the inputs                                                       */
/* ------------------------------------------------------------------------ */

/* The programs that every run needs. */
static const char *const TOOLS[] = { "minimap2", "samtools" };

/* The program that merges a pair of mates, which a run over one file does not need. */
static const char MERGE_TOOL[] = "fastp";

static int check_tools(const mapping_config *cfg, char *error, size_t error_len)
{
    for (size_t i = 0; i < sizeof TOOLS / sizeof *TOOLS; i++) {
        if (!proc_on_path(TOOLS[i])) {
            return fail(TOOLS[i], "is not on PATH", error, error_len);
        }
    }

    if (reads_given(cfg) == MAPPING_MAX_READS && !proc_on_path(MERGE_TOOL)) {
        return fail(MERGE_TOOL, "is not on PATH", error, error_len);
    }

    return 0;
}

/* minimap2 accepts an index in place of a FASTA file. This program needs the bases
 * themselves, because it writes the checksum of every reference into the header, so an
 * index is refused. */
static int check_reference(const mapping_config *cfg, char *error, size_t error_len)
{
    char   magic[sizeof INDEX_MAGIC];
    FILE  *file = fopen(cfg->fasta_path, "rb");
    size_t read;

    if (!file) {
        return fail_with_errno(cfg->fasta_path, error, error_len);
    }

    read = fread(magic, 1, sizeof magic, file);
    fclose(file);

    if (read == sizeof magic && memcmp(magic, INDEX_MAGIC, sizeof magic) == 0) {
        return fail(cfg->fasta_path,
                    "is a minimap2 index; pass the FASTA file it was built from",
                    error, error_len);
    }

    return 0;
}

/* Reads the first bytes of every file of reads to find out what it holds. A format this
 * program cannot align is refused here, because passing it to minimap2 would report the
 * mistake far from where it was made. */
static int classify_reads(const mapping_config *cfg, reads_format *formats,
                          char *error, size_t error_len)
{
    for (size_t i = 0; i < reads_given(cfg); i++) {
        if (access(cfg->reads_paths[i], R_OK) < 0) {
            return fail_with_errno(cfg->reads_paths[i], error, error_len);
        }

        formats[i] = reads_format_of(cfg->reads_paths[i]);

        if (formats[i] == READS_NEITHER) {
            return fail(cfg->reads_paths[i], "is neither a FASTQ nor a BAM file",
                        error, error_len);
        }
    }

    return 0;
}

static int check_pair(const mapping_config *cfg, const reads_format *formats,
                      char *error, size_t error_len)
{
    if (reads_given(cfg) < MAPPING_MAX_READS) {
        return 0;
    }

    if (!preset_accepts_pairs(cfg->preset)) {
        snprintf(error, error_len,
                 "%s reads do not come in pairs; pass one file of reads",
                 preset_name(cfg->preset));
        return -1;
    }

    /* A BAM file holds every read of a run, including both mates of a pair. */
    if (formats[0] == READS_BAM || formats[1] == READS_BAM) {
        return fail_plainly("a BAM file holds both mates; pass it on its own",
                            error, error_len);
    }

    return 0;
}

/* A BAM file that is given on its own must hold raw reads. */
static int check_single_bam(const mapping_config *cfg, const reads_format *formats,
                            char *error, size_t error_len)
{
    const char *path    = cfg->reads_paths[0];
    bool        aligned = false;
    bool        paired  = false;

    if (reads_given(cfg) != 1 || formats[0] != READS_BAM) {
        return 0;
    }

    if (reads_inspect_bam(path, &aligned, &paired, error, error_len) < 0) {
        return -1;
    }

    /* Alignments are what cmuts hmm counts, so a file that holds them belongs to that
     * subcommand and not to this one. */
    if (aligned) {
        return fail(path, "is already aligned; count it with cmuts hmm", error, error_len);
    }

    /* The mates of a pair arrive one after the other, which minimap2 does not read as a
     * pair. */
    if (paired) {
        return fail(path, "holds paired reads; write each mate to a FASTQ file first",
                    error, error_len);
    }

    return 0;
}

/* The sorted alignments are written beside the output file and renamed afterwards, so the
 * directory that will hold the output must already exist. */
static int check_directory(const char *output_path, char *error, size_t error_len)
{
    char        copy[PATH_MAX];
    struct stat info;

    if (snprintf(copy, sizeof copy, "%s", output_path) >= (int)sizeof copy) {
        return fail(output_path, "is too long a path", error, error_len);
    }

    /* The checks run before the pipeline starts its threads.
     * NOLINTNEXTLINE(concurrency-mt-unsafe) */
    const char *directory = dirname(copy);

    if (stat(directory, &info) < 0 || !S_ISDIR(info.st_mode)) {
        return fail(directory, "is not a directory", error, error_len);
    }

    return 0;
}

static int check_output(const mapping_config *cfg, char *error, size_t error_len)
{
    struct stat info;

    if (stat(cfg->output_path, &info) == 0) {
        if (S_ISDIR(info.st_mode)) {
            return fail(cfg->output_path, "is a directory; -o names the file to write",
                        error, error_len);
        }

        if (!cfg->overwrite) {
            return fail(cfg->output_path, "exists; pass --overwrite to replace it",
                        error, error_len);
        }
    }

    return check_directory(cfg->output_path, error, error_len);
}

/* ------------------------------------------------------------------------ */
/* The file that holds the alignments until the run finishes                 */
/* ------------------------------------------------------------------------ */

/* The path is held here as well as by the caller, so that the handler below can remove
 * the file when the run is interrupted. */
static char                  PARTIAL_PATH[PATH_MAX];
static volatile sig_atomic_t PARTIAL_EXISTS;

static const int GUARDED_SIGNALS[] = { SIGINT, SIGTERM, SIGHUP };

#define N_GUARDED_SIGNALS (sizeof GUARDED_SIGNALS / sizeof *GUARDED_SIGNALS)

/* Removes every entry of a directory whose name begins with the prefix. */
static void remove_by_prefix(const char *directory_path, const char *prefix)
{
    /* readdir is safe here, every call being made after the pipeline's threads join.
     * NOLINTNEXTLINE(concurrency-mt-unsafe) */
    DIR   *directory = opendir(directory_path);
    size_t len       = strlen(prefix);
    int    fd;

    if (!directory) {
        return;
    }

    fd = dirfd(directory);

    /* NOLINTBEGIN(concurrency-mt-unsafe) */
    for (const struct dirent *entry = readdir(directory); entry && fd >= 0;
         entry = readdir(directory)) {
        if (strncmp(entry->d_name, prefix, len) == 0) {
            unlinkat(fd, entry->d_name, 0);
        }
    }
    /* NOLINTEND(concurrency-mt-unsafe) */

    closedir(directory);
}

/* Removes the partial output, and every temporary file the sorter wrote beside it.
 *
 * samtools sort writes <partial>.NNNN.bam where the alignments do not fit in memory.
 * Those files remain where it is stopped part way through. Reading a directory is not
 * safe from a signal handler, so only the failure paths call this. */
static void remove_partial_files(const char *partial)
{
    char toward_directory[PATH_MAX];
    char toward_prefix[PATH_MAX];

    unlink(partial);

    if (snprintf(toward_directory, sizeof toward_directory, "%s", partial)
            >= (int)sizeof toward_directory
        || snprintf(toward_prefix, sizeof toward_prefix, "%s", partial)
            >= (int)sizeof toward_prefix) {
        return;
    }

    /* NOLINTNEXTLINE(concurrency-mt-unsafe) */
    remove_by_prefix(dirname(toward_directory), basename(toward_prefix));
}

/* Removes the partial output, then ends this process as the signal would have ended it. */
static void remove_partial(int signal_number)
{
    struct sigaction action = { .sa_handler = SIG_DFL };

    if (PARTIAL_EXISTS) {
        unlink(PARTIAL_PATH);
    }

    sigemptyset(&action.sa_mask);
    sigaction(signal_number, &action, NULL);
    raise(signal_number);
}

/* Names the partial output. The process id separates two runs that write into the same
 * directory. */
static int name_partial(const mapping_config *cfg, char *out, size_t size,
                        char *error, size_t error_len)
{
    if (snprintf(out, size, "%s.partial.%ld", cfg->output_path, (long)getpid())
        >= (int)size) {
        return fail(cfg->output_path, "is too long a path", error, error_len);
    }

    return 0;
}

/* Starts removing the partial output when the run is interrupted.
 *
 * SIGPIPE is ignored as well. A write to a child that has already exited then fails with
 * an error, instead of ending this process before the partial output is removed. */
static void guard_partial(const char *path)
{
    struct sigaction guarded = { .sa_handler = remove_partial };
    struct sigaction ignored = { .sa_handler = SIG_IGN };

    snprintf(PARTIAL_PATH, sizeof PARTIAL_PATH, "%s", path);
    PARTIAL_EXISTS = 1;

    sigemptyset(&guarded.sa_mask);
    sigemptyset(&ignored.sa_mask);

    for (size_t i = 0; i < N_GUARDED_SIGNALS; i++) {
        sigaction(GUARDED_SIGNALS[i], &guarded, NULL);
    }

    sigaction(SIGPIPE, &ignored, NULL);
}

static void unguard_partial(void)
{
    struct sigaction action = { .sa_handler = SIG_DFL };

    PARTIAL_EXISTS = 0;
    sigemptyset(&action.sa_mask);

    for (size_t i = 0; i < N_GUARDED_SIGNALS; i++) {
        sigaction(GUARDED_SIGNALS[i], &action, NULL);
    }

    sigaction(SIGPIPE, &action, NULL);
}

/* ------------------------------------------------------------------------ */
/* The child processes                                                       */
/* ------------------------------------------------------------------------ */

static void write_threads(char *out, size_t size, int threads)
{
    snprintf(out, size, "%d", threads);
}

/* Writes a memory limit as samtools spells one, which is a number and a unit suffix. */
static void write_memory(char *out, size_t size, size_t mebibytes)
{
    snprintf(out, size, "%zuM", mebibytes);
}

/* fastp merges each pair of mates into one read, so that the region where the two mates
 * overlap counts one molecule once. cmuts hmm refuses an unmerged pair for that reason. A
 * pair that cannot be merged is discarded, and the reports are not kept. */
static int spawn_merger(const mapping_config *cfg, int out_fd, pid_t *pid,
                        char *error, size_t error_len)
{
    char threads[THREADS_TEXT_MAX];

    write_threads(threads, sizeof threads,
                  cfg->threads < MERGE_THREADS_MAX ? cfg->threads : MERGE_THREADS_MAX);

    const char *const argv[] = {
        MERGE_TOOL, "--merge", "--stdout",
        "--in1", cfg->reads_paths[0],
        "--in2", cfg->reads_paths[1],
        "--thread", threads,
        "--json", "/dev/null",
        "--html", "/dev/null",
        NULL,
    };

    return proc_spawn(argv, PROC_INHERIT, out_fd, pid, error, error_len);
}

/* input is the file that minimap2 reads, or "-" for its standard input. */
static int spawn_aligner(const mapping_config *cfg, const char *input, int in_fd,
                         int out_fd, pid_t *pid, char *error, size_t error_len)
{
    char threads[THREADS_TEXT_MAX];

    write_threads(threads, sizeof threads, cfg->threads);

    const char *const argv[] = {
        "minimap2", "-a", "-x", preset_name(cfg->preset), "--secondary=no",
        "-t", threads, cfg->fasta_path, input, NULL,
    };

    return proc_spawn(argv, in_fd, out_fd, pid, error, error_len);
}

/* -m is the memory one sorting thread may use, so the total for a run is that figure
 * multiplied by the number of threads. -T is the prefix of the temporary files written
 * each time a thread fills. It is the partial output's own name, so every file one run
 * writes shares a prefix and all of them are removed together. */
static int spawn_sorter(const mapping_config *cfg, const char *partial, int in_fd,
                        pid_t *pid, char *error, size_t error_len)
{
    char threads[THREADS_TEXT_MAX];
    char memory[MEMORY_TEXT_MAX];

    write_threads(threads, sizeof threads, cfg->threads);
    write_memory(memory, sizeof memory, cfg->memory);

    const char *const argv[] = {
        "samtools", "sort", "-@", threads, "-m", memory, "-T", partial,
        "-o", partial, "-", NULL,
    };

    return proc_spawn(argv, in_fd, PROC_INHERIT, pid, error, error_len);
}

/* ------------------------------------------------------------------------ */
/* Correcting the alignments on their way to the sorter                      */
/* ------------------------------------------------------------------------ */

/* Writes the checksum of one reference into the header. A reference that no read was
 * aligned to has no line in the header and is passed over. */
static int write_checksum(sam_hdr_t *header, const cm_fasta_record *record,
                          char *error, size_t error_len)
{
    char hex[CHECKSUM_LEN + 1];

    if (sam_hdr_line_index(header, "SQ", record->name) < 0) {
        return 0;
    }

    if (!checksum_sequence(record->seq, record->len, hex)) {
        snprintf(error, error_len, "the checksum of reference \"%s\" could not be "
                                   "computed", record->name);
        return -1;
    }

    if (sam_hdr_update_line(header, "SQ", "SN", record->name, "M5", hex, NULL) < 0) {
        snprintf(error, error_len, "the checksum of reference \"%s\" could not be "
                                   "written to the header", record->name);
        return -1;
    }

    return 0;
}

/* Writes the checksum of every reference into the header.
 *
 * minimap2 writes the name and the length of each reference. cmuts hmm compares the
 * checksum in the header against the FASTA file it is given, which confirms that the
 * FASTA file holds the sequences the reads were aligned to. */
static int write_checksums(sam_hdr_t *header, const char *fasta_path,
                           char *error, size_t error_len)
{
    const char      *why    = NULL;
    cm_fasta_reader *reader = cm_fasta_open(fasta_path, &why);
    cm_fasta_record  record;
    int              status = 0;

    if (!reader) {
        return fail(fasta_path, why, error, error_len);
    }

    while (status == 0 && cm_fasta_next(reader, &record) == CM_ITER_OK) {
        status = write_checksum(header, &record, error, error_len);
    }

    if (status == 0 && cm_fasta_error(reader)) {
        status = fail(fasta_path, cm_fasta_error(reader), error, error_len);
    }

    cm_fasta_close(reader);

    return status;
}

/* Writes one alignment out, shortening any read name that would not fit in a BAM file.
 *
 * minimap2 copies the name from the reads, and neither a FASTQ file nor a BAM file limits
 * its length. The alignments are handled as text because htslib refuses to parse a record
 * whose name is already too long. */
static int write_record_line(FILE *out, const char *line, size_t len)
{
    const char *tab  = memchr(line, '\t', len);
    size_t      name = tab ? (size_t)(tab - line) : len;
    size_t      rest = len - name;
    size_t      keep = name < NAME_LIMIT ? name : NAME_LIMIT;

    if (fwrite(line, 1, keep, out) != keep) {
        return -1;
    }

    return fwrite(line + name, 1, rest, out) == rest ? 0 : -1;
}

/* The alignments arrive from minimap2 and leave for samtools sort. */
typedef struct {
    FILE      *in;
    FILE      *out;
    char      *line;
    size_t     line_size;
    ssize_t    pending;  /* the length of the first line after the header */
    sam_hdr_t *header;
} corrector;

/* Reads every line at the start of the stream that begins with an at sign. The first line
 * after the header stays in the line buffer. */
static int read_header(corrector *c, char *error, size_t error_len)
{
    kstring_t text = KS_INITIALIZE;
    ssize_t   len;

    while ((len = getline(&c->line, &c->line_size, c->in)) > 0) {
        if (c->line[0] != '@') {
            c->pending = len;
            break;
        }

        kputsn(c->line, (size_t)len, &text);
    }

    if (text.l == 0) {
        ks_free(&text);
        return fail("minimap2", "wrote no header", error, error_len);
    }

    c->header = sam_hdr_parse(text.l, text.s);
    ks_free(&text);

    if (!c->header) {
        return fail("minimap2", "wrote a header that could not be read", error, error_len);
    }

    return 0;
}

/* Writes the header on to the sorter, with the checksum of every reference added. */
static int write_header(corrector *c, const mapping_config *cfg,
                        char *error, size_t error_len)
{
    if (write_checksums(c->header, cfg->fasta_path, error, error_len) < 0) {
        return -1;
    }

    if (fputs(sam_hdr_str(c->header), c->out) == EOF) {
        return fail("samtools sort", "did not accept the header", error, error_len);
    }

    return 0;
}

static int correct_records(corrector *c, char *error, size_t error_len)
{
    ssize_t len = c->pending;

    while (len > 0) {
        if (write_record_line(c->out, c->line, (size_t)len) < 0) {
            return fail("samtools sort", "did not accept an alignment", error, error_len);
        }

        len = getline(&c->line, &c->line_size, c->in);
    }

    if (ferror(c->in)) {
        return fail("minimap2", "could not be read", error, error_len);
    }

    return 0;
}

/* Corrects everything that arrives on one stream and writes it to the other. */
static int correct_stream(const mapping_config *cfg, FILE *in, FILE *out,
                          char *error, size_t error_len)
{
    corrector c      = { .in = in, .out = out, .line_size = LINE_SIZE };
    int       status;

    c.line = malloc(c.line_size);

    if (!c.line) {
        return fail_plainly("out of memory", error, error_len);
    }

    status = read_header(&c, error, error_len);

    if (status == 0) {
        status = write_header(&c, cfg, error, error_len);
    }

    if (status == 0) {
        /* getline grows the line buffer, which the analyzer reads as a second allocation
         * that the free below does not reach.
         * NOLINTNEXTLINE(clang-analyzer-unix.Malloc) */
        status = correct_records(&c, error, error_len);
    }

    free(c.line);
    sam_hdr_destroy(c.header);

    return status;
}

/* Reads the alignments from the aligner and passes them to the sorter, correcting the
 * header and the read names on the way. Both descriptors are closed here. */
static int correct_alignments(const mapping_config *cfg, int in_fd, int out_fd,
                              char *error, size_t error_len)
{
    FILE *in  = fdopen(in_fd, "r");
    FILE *out = in ? fdopen(out_fd, "w") : NULL;
    int   status;

    if (!in || !out) {
        close_either(in, in_fd);
        close_either(out, out_fd);

        return fail_plainly("the aligner and the sorter could not be connected",
                            error, error_len);
    }

    status = correct_stream(cfg, in, out, error, error_len);
    fclose(in);

    /* The last of the alignments reach the sorter as this stream is flushed. */
    if (fclose(out) != 0 && status == 0) {
        status = fail("samtools sort", "did not accept the alignments", error, error_len);
    }

    return status;
}

/* ------------------------------------------------------------------------ */
/* Unpacking a BAM file                                                      */
/* ------------------------------------------------------------------------ */

/* A BAM file is written out as FASTQ on a thread of its own. This process reads the
 * output of the aligner at the same time, and one thread doing both would stop as soon as
 * either pipe filled up. */
typedef struct {
    const char *path;
    int         fd;
    int         status;
    char        error[CM_ERROR_MAX];
} unpacker;

static void *unpack(void *given)
{
    unpacker *u = given;

    u->status = reads_write_fastq(u->path, u->fd, u->error, sizeof u->error);

    return NULL;
}

/* ------------------------------------------------------------------------ */
/* Running the pipeline                                                      */
/* ------------------------------------------------------------------------ */

/* Every child process of one run, and the thread that supplies the aligner's input. */
typedef struct {
    pid_t     merger;  /* fastp, or 0 when the reads are not a pair */
    pid_t     aligner;
    pid_t     sorter;
    pthread_t thread;
    bool      unpacking;
    unpacker  work;
} stages;

static void close_if_open(int *fd)
{
    if (*fd >= 0) {
        close(*fd);
        *fd = -1;
    }
}

/* Starts the stage that supplies the aligner's input, and returns the file the aligner
 * reads.
 *
 * A pair of files is merged by fastp. A BAM file is unpacked by the thread above. A
 * single FASTQ file is passed to the aligner by name, so no stage is started and its
 * standard input is unchanged. Nothing is left open where this fails. */
static int start_source(const mapping_config *cfg, const reads_format *formats,
                        stages *s, const char **input, int *aligner_stdin,
                        char *error, size_t error_len)
{
    int fds[2];

    *input         = cfg->reads_paths[0];
    *aligner_stdin = PROC_INHERIT;

    if (reads_given(cfg) == 1 && formats[0] == READS_FASTQ) {
        return 0;
    }

    if (proc_pipe(fds, error, error_len) < 0) {
        return -1;
    }

    *input         = "-";
    *aligner_stdin = fds[0];

    if (reads_given(cfg) == MAPPING_MAX_READS) {
        if (spawn_merger(cfg, fds[1], &s->merger, error, error_len) < 0) {
            close_if_open(aligner_stdin);
            return -1;
        }

        return 0;
    }

    s->work.path = cfg->reads_paths[0];
    s->work.fd   = fds[1];

    if (pthread_create(&s->thread, NULL, unpack, &s->work) != 0) {
        close(fds[1]);
        close_if_open(aligner_stdin);

        return fail(cfg->reads_paths[0], "could not be unpacked", error, error_len);
    }

    s->unpacking = true;

    return 0;
}

/* Starts every child, and returns the two descriptors that this process keeps: the one
 * the alignments arrive on, and the one they leave by. */
static int start_stages(const mapping_config *cfg, const reads_format *formats,
                        const char *partial, stages *s, int *from_aligner, int *to_sorter,
                        char *error, size_t error_len)
{
    const char *input         = NULL;
    int         aligner_stdin = PROC_INHERIT;
    int         aligned[2];
    int         sorted[2];

    if (start_source(cfg, formats, s, &input, &aligner_stdin, error, error_len) < 0) {
        return -1;
    }

    if (proc_pipe(aligned, error, error_len) < 0) {
        close_if_open(&aligner_stdin);
        return -1;
    }

    *from_aligner = aligned[0];

    /* The spawn closes the descriptors it passes on, whether or not it succeeds. */
    if (spawn_aligner(cfg, input, aligner_stdin, aligned[1], &s->aligner,
                      error, error_len) < 0) {
        return -1;
    }

    if (proc_pipe(sorted, error, error_len) < 0) {
        return -1;
    }

    *to_sorter = sorted[1];

    return spawn_sorter(cfg, partial, sorted[0], &s->sorter, error, error_len);
}

/* Records a failure if nothing has failed yet, and reports a failure either way. */
static int keep_first(int status, const char *why, char *error, size_t error_len)
{
    if (status == 0) {
        snprintf(error, error_len, "%s", why);
    }

    return -1;
}

/* Waits for the thread and for every child that was started. The first failure is the one
 * that is reported, because the others usually follow from it. */
static int finish_stages(stages *s, int status, char *error, size_t error_len)
{
    char scratch[CM_ERROR_MAX];

    if (s->unpacking) {
        pthread_join(s->thread, NULL);

        if (s->work.status < 0) {
            status = keep_first(status, s->work.error, error, error_len);
        }
    }

    if (s->merger != 0 && proc_wait(s->merger, MERGE_TOOL, scratch, sizeof scratch) < 0) {
        status = keep_first(status, scratch, error, error_len);
    }

    if (s->aligner != 0
        && proc_wait(s->aligner, "minimap2", scratch, sizeof scratch) < 0) {
        status = keep_first(status, scratch, error, error_len);
    }

    if (s->sorter != 0
        && proc_wait(s->sorter, "samtools sort", scratch, sizeof scratch) < 0) {
        status = keep_first(status, scratch, error, error_len);
    }

    return status;
}

static int run_stages(const mapping_config *cfg, const reads_format *formats,
                      const char *partial, char *error, size_t error_len)
{
    stages s            = { 0 };
    int    from_aligner = -1;
    int    to_sorter    = -1;
    int    status       = start_stages(cfg, formats, partial, &s, &from_aligner,
                                       &to_sorter, error, error_len);

    if (status == 0) {
        status       = correct_alignments(cfg, from_aligner, to_sorter, error, error_len);
        from_aligner = -1;
        to_sorter    = -1;
    }

    /* Closing these ends the input of a child that is still running. */
    close_if_open(&from_aligner);
    close_if_open(&to_sorter);

    return finish_stages(&s, status, error, error_len);
}

/* ------------------------------------------------------------------------ */
/* The run                                                                   */
/* ------------------------------------------------------------------------ */

int mapping_run(const mapping_config *cfg, char *error, size_t error_len)
{
    reads_format formats[MAPPING_MAX_READS];
    char         partial[PATH_MAX];
    int          status;

    if (check_tools(cfg, error, error_len) < 0
        || check_reference(cfg, error, error_len) < 0
        || classify_reads(cfg, formats, error, error_len) < 0
        || check_pair(cfg, formats, error, error_len) < 0
        || check_single_bam(cfg, formats, error, error_len) < 0
        || check_output(cfg, error, error_len) < 0
        || name_partial(cfg, partial, sizeof partial, error, error_len) < 0) {
        return -1;
    }

    guard_partial(partial);

    status = run_stages(cfg, formats, partial, error, error_len);

    if (status == 0 && rename(partial, cfg->output_path) < 0) {
        status = fail_with_errno(cfg->output_path, error, error_len);
    }

    if (status < 0) {
        remove_partial_files(partial);
    }

    unguard_partial();

    return status;
}
