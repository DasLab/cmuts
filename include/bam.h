/* bam.h -- sequential iteration over the alignments of a BAM file.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#pragma once

#include <stdbool.h>
#include <stdint.h>

#include <htslib/sam.h>

#include "iter.h"

/* A single alignment record.
 *
 * Every pointer borrows memory owned by the reader that produced the record,
 * and is invalidated by the next cm_bam_next() or cm_bam_close() call on that
 * reader. Copy anything that must outlive the current iteration step.
 *
 * No filtering is applied, so unmapped reads reach the caller: those carry
 * tid == CM_BAM_NO_REF, pos == -1 and n_cigar == 0. Test flag against
 * BAM_FUNMAP rather than tid, since a paired read may be unmapped and still
 * carry its mate's reference. */
typedef struct {
    int32_t         tid;      /* reference index into the file header */
    hts_pos_t       pos;      /* 0-based leftmost coordinate on that reference */
    uint16_t        flag;     /* SAM flags, for BAM_FUNMAP and friends */
    uint8_t         mapq;     /* mapping quality */
    const uint32_t *cigar;    /* raw ops; bam_cigar_op / bam_cigar_oplen decode */
    uint32_t        n_cigar;  /* number of ops in cigar */
    const char     *md;       /* MD tag; NULL when the aligner omitted it */
    const uint8_t  *qual;     /* raw PHRED scores, no +33 offset; NULL if absent */
    int32_t         l_qseq;   /* number of bases, and length of qual */
} cm_bam_record;

/* Reference index of an unmapped read, and of any tid with no header entry. */
#define CM_BAM_NO_REF (-1)

typedef struct cm_bam_reader cm_bam_reader;

/* Opens path and reads its header. No .bai index is required or consulted:
 * iteration is a linear pass, in file order. Returns NULL on failure, with
 * errno set by the underlying open. */
cm_bam_reader *cm_bam_open(const char *path);

/* Fills out with the next alignment. Returns a cm_iter_status; on
 * CM_ITER_ERROR the cause is available from cm_bam_error(). */
int cm_bam_next(cm_bam_reader *reader, cm_bam_record *out);

void cm_bam_close(cm_bam_reader *reader);

/* Description of the reader's failure, or NULL if it has not failed. */
const char *cm_bam_error(const cm_bam_reader *reader);

/* Adds worker threads for BGZF decompression, which parallelizes inflation
 * only: reading and record parsing stay sequential. Worth setting when the
 * loader is the bottleneck, and pointless on small files. */
int cm_bam_set_threads(cm_bam_reader *reader, int threads);

/* ------------------------------------------------------------------------ */
/* The record underlying a view                                              */
/* ------------------------------------------------------------------------ */

/* The htslib record behind the most recent cm_bam_next(), for callers that
 * need to retain a read past the next advance. Copy it with bam_copy1. */
const bam1_t *cm_bam_raw(const cm_bam_reader *reader);

/* Presents an htslib record as a cm_bam_record. The view borrows from record,
 * so it stays valid exactly as long as record does. */
void cm_bam_record_view(const bam1_t *record, cm_bam_record *out);

/* ------------------------------------------------------------------------ */
/* Header queries                                                            */
/* ------------------------------------------------------------------------ */

/* Name of the reference with the given index, or NULL if there is no such
 * reference. The returned string is owned by the reader's header. */
const char *cm_bam_refname(const cm_bam_reader *reader, int32_t tid);

/* Length of the reference with the given index, or -1 if there is no such
 * reference. */
hts_pos_t cm_bam_reflen(const cm_bam_reader *reader, int32_t tid);

/* Length of the longest reference in the header. */
hts_pos_t cm_bam_max_reflen(const cm_bam_reader *reader);

/* Number of references declared in the header. */
int32_t cm_bam_nref(const cm_bam_reader *reader);

/* Whether the header declares SO:coordinate. Reference-grouped input is what
 * lets a reference be finished the moment the reader moves past it; without
 * it, every reference stays live to the end of the file. */
bool cm_bam_is_coordinate_sorted(const cm_bam_reader *reader);
