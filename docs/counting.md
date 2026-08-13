# How cmuts counts

## Reads that never count

Unmapped reads are counted separately and reach no reference. Secondary alignments, records storing no sequence, and records carrying no CIGAR are always refused, as is a read of mapping quality 255, which the SAM specification defines as "unavailable" rather than as a score above every threshold.

Every mapped read is either counted or rejected, and both totals are written, so nothing goes missing between the file and the result.

## What counts as a modification

A substitution, a deletion and an insertion each carry a weight between 0 and 1, set by `--substitution-weight`, `--deletion-weight` and `--insertion-weight`. A weight of zero leaves that kind of difference out of the total entirely.

Insertions are weighed at zero by default. An inserted base sits between two reference positions rather than at one, so counting it means choosing which of them it belongs to.

## Gaps that could be written anywhere

A deletion inside a run of the same base can be written after any base of the run, and every CIGAR that writes it describes the same alignment. cmuts marginalizes over the placements the aligner might have chosen instead, so a result does not depend on which one it did choose.

`--band` is how far either side of the CIGAR the marginal looks, in reference positions. It must be at least as wide as a gap for the placements of that gap to agree. The default of 2 covers the gaps that occur in practice, and a wider band costs time.

## Positions that get a rate

`--min-depth` is the evidence a position needs before a rate is written for it. Below it the reactivity and its error are missing rather than zero: a position nothing was observed at has no rate, which says something different from a position where nothing was modified.

The default of 1 is one whole observation. Below that, the standard error of a proportion is divided by a fraction and stops being bounded by a half, so a rate reported there would carry an error that says little.

## What subtraction does to each dataset

<!-- BEGIN GENERATED RULES -->
| Dataset | Without a control | With one |
| --- | --- | --- |
| `coverage` | every input added | the same |
| `reactivity` | the untreated value taken from the treated one | that difference over the control's value, and missing where the control measured nothing to divide by |
| `error` | the two added in quadrature | the error of that ratio, taking the three runs as independent |
| `reads/lengths` | every input added | the same |
| `reads/counted` | every input added | the same |
| `reads/rejected` | every input added | the same |
| `reads/unmapped` | every input added | the same |
<!-- END GENERATED RULES -->

A rate is missing wherever either input is missing one, since a difference needs both. `--clip` raises a negative difference to zero and leaves a missing value missing, rather than reading it as a zero to raise.
