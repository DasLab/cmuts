/* simulate.c -- alignments built from a reference and a mutation model.
 *
 * A read is first laid out as a list of single-base events, from which the CIGAR, SEQ, MD and NM
 * are all then derived. Nothing is written twice, so the four cannot contradict each other
 * however the events fall.
 *
 * Author: Hamish M. Blair <hmblair@stanford.edu>
 */

#include "simulate.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char BASES[] = "ACGT";

/* How much room an MD string is given: a few characters per event, plus room for the
 * count it closes with.
 * Neither figure is what keeps the writing inside the buffer -- the cursor is held
 * short of the end at every step regardless -- so an underestimate costs a tag cut
 * short and nothing worse. MD_HEADROOM is the margin the loop keeps so that a tag cut
 * short ends after a whole event and not partway through a number. */
#define MD_PER_EVENT 4
#define MD_TAIL      32
#define MD_HEADROOM  16

typedef enum {
    EV_MATCH,
    EV_MISMATCH,
    EV_INSERT,
    EV_DELETE,
    EV_SOFT,
} ev_kind;

typedef struct {
    unsigned char kind;
    char          ref;    /* reference base, for MATCH, MISMATCH and DELETE */
    char          query;  /* read base, for everything except DELETE */
} sim_event;

struct sim_scratch {
    sim_event *events;
    size_t     capacity;
    uint32_t  *cigar;
    char      *seq;
    char      *qual;
    char      *md;
    size_t     md_capacity;
};

/* ------------------------------------------------------------------------ */
/* Scratch                                                                   */
/* ------------------------------------------------------------------------ */

static size_t largest_read(const sim_model *model)
{
    size_t span   = (size_t)distribution_maximum(&model->length);
    size_t insert = (size_t)distribution_maximum(&model->insertions) *
                    (size_t)distribution_maximum(&model->insertion_length);
    size_t clip   = 2 * (size_t)distribution_maximum(&model->soft_clip_length);

    return span + insert + clip + 1;
}

sim_scratch *sim_scratch_create(const sim_model *model)
{
    sim_scratch *scratch = calloc(1, sizeof *scratch);
    if (!scratch) {
        return NULL;
    }

    scratch->capacity    = largest_read(model);
    scratch->md_capacity = scratch->capacity * MD_PER_EVENT + MD_TAIL;

    scratch->events = calloc(scratch->capacity, sizeof *scratch->events);
    scratch->cigar  = calloc(scratch->capacity, sizeof *scratch->cigar);
    scratch->seq    = calloc(scratch->capacity + 1, 1);
    scratch->qual   = calloc(scratch->capacity + 1, 1);
    scratch->md     = calloc(scratch->md_capacity, 1);

    if (!scratch->events || !scratch->cigar || !scratch->seq ||
        !scratch->qual || !scratch->md) {
        sim_scratch_destroy(scratch);
        return NULL;
    }

    return scratch;
}

void sim_scratch_destroy(sim_scratch *scratch)
{
    if (!scratch) {
        return;
    }

    free(scratch->md);
    free(scratch->qual);
    free(scratch->seq);
    free(scratch->cigar);
    free(scratch->events);
    free(scratch);
}

/* ------------------------------------------------------------------------ */
/* Bases                                                                     */
/* ------------------------------------------------------------------------ */

static char random_base(rng *r)
{
    return BASES[rng_next(r) & 3];
}

static char different_base(char from, rng *r)
{
    char base;

    do {
        base = random_base(r);
    } while (base == from);

    return base;
}

void sim_sequence(char *out, size_t len, rng *r)
{
    for (size_t i = 0; i < len; i++) {
        out[i] = random_base(r);
    }

    out[len] = '\0';
}

/* ------------------------------------------------------------------------ */
/* Laying out a read                                                         */
/* ------------------------------------------------------------------------ */

static size_t add_clip(sim_event *events, size_t n, size_t cap, long len, rng *r)
{
    for (long i = 0; i < len && n < cap; i++) {
        events[n].kind  = EV_SOFT;
        events[n].query = random_base(r);
        n++;
    }

    return n;
}

/* Returns whether an event falls at this position, spreading the remaining events
 * uniformly over the positions left so that no list of positions has to be drawn and
 * sorted first. */
static bool event_falls_here(long remaining, size_t positions_left, rng *r)
{
    return remaining > 0 &&
           rng_chance(r, (double)remaining / (double)positions_left);
}

sim_placement sim_place(const sim_model *model, rng *r, size_t reflen)
{
    sim_placement where;
    long          span = distribution_draw(&model->length, r);

    if (span < 1) {
        span = 1;
    }
    if ((size_t)span > reflen) {
        span = (long)reflen;
    }

    where.span  = (size_t)span;
    where.start = rng_between(r, 0, (long)(reflen - where.span));

    return where;
}

static size_t lay_out_read(sim_scratch *scratch, const sim_model *model, rng *r,
                           const char *reference, sim_placement out)
{
    sim_event *events = scratch->events;
    size_t     cap    = scratch->capacity;

    long clips     = distribution_draw(&model->soft_clips, r);
    long insertion = distribution_draw(&model->insertions, r);
    long deletion  = distribution_draw(&model->deletions, r);
    size_t n       = 0;

    if (clips > 0) {
        n = add_clip(events, n, cap,
                     distribution_draw(&model->soft_clip_length, r), r);
    }

    for (size_t i = 0; i < out.span && n < cap; ) {
        size_t left = out.span - i;

        /* Indels are kept away from either end, so that the CIGAR neither begins nor ends with
         * an operation consuming only one of the two sequences. */
        bool interior = i > 0 && left > 1;

        if (interior && event_falls_here(deletion, left, r)) {
            long len = distribution_draw(&model->deletion_length, r);

            if (len > (long)left - 1) {
                len = (long)left - 1;
            }

            if (len > 0) {
                for (long k = 0; k < len && n < cap; k++) {
                    events[n].kind = EV_DELETE;
                    events[n].ref  = reference[out.start + (long)i + k];
                    n++;
                }
                i += (size_t)len;
                deletion--;
                continue;
            }
        }

        if (interior && event_falls_here(insertion, left, r)) {
            long len = distribution_draw(&model->insertion_length, r);

            for (long k = 0; k < len && n < cap; k++) {
                events[n].kind  = EV_INSERT;
                events[n].query = random_base(r);
                n++;
            }
            insertion--;
        }

        if (n >= cap) {
            break;
        }

        char base     = reference[out.start + i];
        bool mismatch = rng_chance(r, model->mismatch_rate);

        events[n].kind  = mismatch ? EV_MISMATCH : EV_MATCH;
        events[n].ref   = base;
        events[n].query = mismatch ? different_base(base, r) : base;
        n++;
        i++;
    }

    if (clips > 1) {
        n = add_clip(events, n, cap,
                     distribution_draw(&model->soft_clip_length, r), r);
    }

    return n;
}

/* ------------------------------------------------------------------------ */
/* Deriving the encodings                                                    */
/* ------------------------------------------------------------------------ */

static int cigar_op_of(unsigned char kind)
{
    switch (kind) {
        case EV_MATCH:
        case EV_MISMATCH: return BAM_CMATCH;
        case EV_INSERT:   return BAM_CINS;
        case EV_DELETE:   return BAM_CDEL;
        default:          return BAM_CSOFT_CLIP;
    }
}

static size_t build_cigar(const sim_event *events, size_t n, uint32_t *out)
{
    size_t written = 0;

    for (size_t i = 0; i < n; ) {
        int    op  = cigar_op_of(events[i].kind);
        size_t run = 1;

        while (i + run < n && cigar_op_of(events[i + run].kind) == op) {
            run++;
        }

        out[written++] = bam_cigar_gen((uint32_t)run, op);
        i += run;
    }

    return written;
}

static size_t build_sequence(const sim_event *events, size_t n, const sim_model *model,
                             rng *r, char *seq, char *qual)
{
    size_t len = 0;

    for (size_t i = 0; i < n; i++) {
        if (events[i].kind == EV_DELETE) {
            continue;
        }

        seq[len]  = events[i].query;
        qual[len] = (char)distribution_draw(&model->base_quality, r);
        len++;
    }

    seq[len] = '\0';
    return len;
}

/* Appends a match count and the one character closing it: the base that broke the run, or the
 * caret opening a deletion.
 *
 * Returns the cursor where the text ends, not where it would have ended had there been room,
 * which is what snprintf reports. The two differ only once the buffer is full, and a cursor past
 * the end would put every write after it outside the buffer. */
static size_t append_run(char *md, size_t used, size_t cap, long run,
                         char closing)
{
    int wrote = snprintf(md + used, cap - used, "%ld%c", run, closing);

    if (wrote < 0) {
        return used;
    }

    return used + (size_t)wrote < cap ? used + (size_t)wrote : cap - 1;
}

/* Appends one deleted base, written only where both it and the terminator after it fit. */
static size_t append_base(char *md, size_t used, size_t cap, char base)
{
    if (used + 1 >= cap) {
        return used;
    }

    md[used++] = base;
    md[used]   = '\0';

    return used;
}

/* Builds the MD tag. It describes only the reference-consuming positions, so insertions and soft
 * clips contribute nothing. It opens and closes with a match count and carries one between every
 * pair of events, which is why a run of zero is still written out.
 *
 * Every step leaves the cursor short of cap, so the closing count writes inside the buffer
 * however the events fell. */
static void build_md(const sim_event *events, size_t n, char *md, size_t cap)
{
    size_t used = 0;
    long   run  = 0;
    bool   deleting = false;

    for (size_t i = 0; i < n && used + MD_HEADROOM < cap; i++) {
        switch (events[i].kind) {
            case EV_MATCH:
                deleting = false;
                run++;
                break;

            case EV_MISMATCH:
                deleting = false;
                used = append_run(md, used, cap, run, events[i].ref);
                run  = 0;
                break;

            case EV_DELETE:
                if (!deleting) {
                    used = append_run(md, used, cap, run, '^');
                    run  = 0;
                    deleting = true;
                }
                used = append_base(md, used, cap, events[i].ref);
                break;

            default:
                break;
        }
    }

    snprintf(md + used, cap - used, "%ld", run);
}

static int count_differences(const sim_event *events, size_t n)
{
    int nm = 0;

    for (size_t i = 0; i < n; i++) {
        if (events[i].kind != EV_MATCH && events[i].kind != EV_SOFT) {
            nm++;
        }
    }

    return nm;
}

/* ------------------------------------------------------------------------ */
/* Records                                                                   */
/* ------------------------------------------------------------------------ */

static int attach_tags(bam1_t *rec, int nm, const char *md)
{
    int32_t value = nm;

    if (bam_aux_append(rec, "NM", 'i', sizeof value, (const uint8_t *)&value) < 0) {
        return -1;
    }

    return bam_aux_append(rec, "MD", 'Z', (int)strlen(md) + 1,
                          (const uint8_t *)md) < 0 ? -1 : 0;
}

int sim_alignment(bam1_t *rec, sim_scratch *scratch, const sim_model *model,
                  rng *r, const char *name, int32_t tid,
                  const char *reference, sim_placement where)
{
    size_t   events  = lay_out_read(scratch, model, r, reference, where);
    size_t   n_cigar = build_cigar(scratch->events, events, scratch->cigar);
    size_t   l_seq   = build_sequence(scratch->events, events, model, r,
                                      scratch->seq, scratch->qual);
    uint16_t flag    = rng_chance(r, model->reverse_fraction) ? BAM_FREVERSE : 0;
    uint8_t  mapq    = (uint8_t)distribution_draw(&model->mapq, r);

    build_md(scratch->events, events, scratch->md, scratch->md_capacity);

    if (bam_set1(rec, strlen(name), name, flag, tid, where.start, mapq,
                 n_cigar, scratch->cigar, -1, -1, 0,
                 l_seq, scratch->seq, scratch->qual, 0) < 0) {
        return -1;
    }

    return attach_tags(rec, count_differences(scratch->events, events), scratch->md);
}

int sim_unmapped(bam1_t *rec, sim_scratch *scratch, const sim_model *model,
                 rng *r, const char *name)
{
    long len = distribution_draw(&model->length, r);

    if (len < 1) {
        len = 1;
    }
    if ((size_t)len > scratch->capacity) {
        len = (long)scratch->capacity;
    }

    sim_sequence(scratch->seq, (size_t)len, r);
    for (long i = 0; i < len; i++) {
        scratch->qual[i] = (char)distribution_draw(&model->base_quality, r);
    }

    return bam_set1(rec, strlen(name), name, BAM_FUNMAP, -1, -1, 0,
                    0, NULL, -1, -1, 0,
                    (size_t)len, scratch->seq, scratch->qual, 0) < 0 ? -1 : 0;
}
