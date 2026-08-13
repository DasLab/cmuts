# Subtracting a background

An untreated sample carries the mutations that were there without the probe.
`cmuts-sub` takes one cmuts output off another:

```sh
cmuts-sub -o difference.h5 treated.h5 untreated.h5
```

Both inputs must describe the same references, at the same width. The result
holds the same datasets as its inputs, so anything that reads a cmuts output
reads it too.

## What happens to each dataset

| Dataset | Result |
| --- | --- |
| `reactivity` | the treated rate less the untreated one |
| `error` | the two errors added in quadrature |
| `coverage` | the two coverages added |
| `reads/*` | the two counts added |

A rate is missing wherever either input is missing one, since a difference
needs both.

## A negative reactivity

Where the background is above the signal, the difference is negative. That is
the measurement, and it is kept. `--clip` raises a negative difference to zero
instead, for the downstream tools that require it.

Clipping does not touch a missing value: a position with no rate keeps none,
rather than becoming zero.

## A denatured control

A denatured sample is one where every position should be reactive, so it
measures how reachable each position is regardless of structure. Dividing by it
puts profiles from different references on a comparable scale:

```sh
cmuts-sub -o normalized.h5 -d denatured.h5 treated.h5 untreated.h5
```

The reactivity becomes the background-subtracted rate over the control's, and
the error follows the division. Where the control measured nothing at a
position, there is nothing to divide by and the result is missing there.
