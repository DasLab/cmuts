# Running cmuts

```sh
cmuts -f references.fasta -o results.h5 alignments.bam
```

The alignments must be coordinate sorted, and cmuts refuses a file that is not. `samtools sort` produces one. The FASTA must hold the same references, in the same order as the alignment header declares them; where the header carries an M5 checksum, cmuts checks the sequences against it and refuses a FASTA that does not match.

An existing output file is never replaced without `--overwrite`.

## Which reads are counted

Unmapped reads are counted separately and never reach a reference. Secondary alignments, records storing no sequence, and records carrying no CIGAR are always refused. So is a read of mapping quality 255, which the SAM specification defines as "unavailable" rather than as a score above every threshold.

Beyond that, `--min-mapq`, `--min-length`, `--max-length` and `--strand` decide what is kept. The length bounds are on the sequence a read stores, which is not the span it covers on the reference: a read carrying insertions or soft-clipped ends stores more than it aligns to, and one carrying deletions stores less.

Every mapped read is either counted or rejected, and both totals are written, so nothing goes missing between the file and the result.

## What is counted as a modification

A substitution, a deletion and an insertion each carry a weight between 0 and 1, set by `--substitution-weight`, `--deletion-weight` and `--insertion-weight`. A weight of zero leaves that kind of difference out of the total entirely.

Insertions are weighed at zero by default. An inserted base sits between two reference positions rather than at one, so counting it means deciding which position it belongs to.

## Ambiguous gaps

A deletion inside a run of the same base can be written after any base of the run, and every CIGAR that writes it describes the same alignment. cmuts marginalizes over the placements the aligner might have chosen instead, so the result does not depend on which one it did choose.

`--band` is how far either side of the CIGAR the marginal looks, in reference positions. It must be at least as wide as the gap for the placements of that gap to agree; the default of 2 covers the gaps that occur in practice, and a wider band costs time.

## Which positions get a rate

`--min-depth` is the evidence a position needs before a rate is written for it. Below it, the reactivity and its error are missing rather than zero: a position nothing was observed at has no rate, which is a different statement from a position where nothing was modified.

The default of 1 is one whole observation. Below that the standard error of a proportion is divided by a fraction and stops being bounded by a half, so a rate reported there would carry an error that says little.
