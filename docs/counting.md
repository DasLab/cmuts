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

| Dataset | After `cmuts-sub` |
| --- | --- |
| `reactivity` | the treated rate less the untreated one |
| `error` | the two errors added in quadrature |
| `coverage` | the two coverages added |
| `reads/*` | the two counts added |

A rate is missing wherever either input is missing one, since a difference needs both. `--clip` raises a negative difference to zero and leaves a missing value missing, rather than reading it as a zero to raise.

## What division does to each dataset

| Dataset | After `cmuts-div` |
| --- | --- |
| `reactivity` | the rate over the control's |
| `error` | the two errors added in quadrature, the control's scaled by the rate, over the control's rate |
| `coverage` | the two coverages added |
| `reads/*` | the two counts added |

Where the control measured nothing at a position, there is nothing to divide by and the rate and its error are missing there.

## What normalization does to each dataset

| Dataset | After `cmuts-norm` |
| --- | --- |
| `reactivity` | the rate over the scale, then held within the clipping bounds |
| `error` | the error over the same scale |
| `coverage` | unchanged |
| `reads/*` | unchanged |

The scale is one number drawn from the rates of every input given to the run, so a rate is comparable across those inputs and not across separate runs. A scale that comes out as zero, negative, or undefined leaves the rates as they are.
