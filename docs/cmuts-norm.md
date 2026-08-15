# cmuts-norm

## Purpose

Normalization of reactivity rates against a scale taken from the rates themselves.

## Requires

- One or more sets of reactivity rates, in cmuts-compatible HDF5 files

## Usage

Every input given to one run shares a single scale, and each is written to an output of its own. `--output` is repeated once per input and paired with them in order.

```sh
cmuts-norm -o apo.norm.h5 -o holo.norm.h5 apo.h5 holo.h5
```

Running the inputs separately gives each its own scale instead.

```sh
cmuts-norm -o apo.norm.h5 apo.h5
cmuts-norm -o holo.norm.h5 holo.h5
```

Pool the inputs whenever the rates are to be compared across experiments, since a scale of its own puts each experiment on a different footing. There is no ordering requirement against [cmuts-sub](cmuts-sub.md): normalizing before subtraction and normalizing after it differ only in which rates the scale is drawn from.

## The Scale

The scale is one number, so the inputs need not have been counted against the same references.

`--norm ubr` (the default) takes the 90th percentile of the pooled rates, counting only positions whose coverage exceeds `--min-coverage`. `--norm outlier` drops the highest 2% of the pooled rates as outliers and averages what lies between there and the highest 10%; it reads no coverage, so `--min-coverage` does not apply to it.

Both divide the reactivity and its error by the scale, leaving the coverage and the four read-related datasets alone. A scale that comes out as zero, negative, or undefined leaves the rates as they are. Each output records the scale it was given as a `norm` attribute on its root group.

```{note}
`--min-coverage` defaults to 500, which was written against a single library's coverage. A [cmuts-sub](cmuts-sub.md) output holds the sum of the treated and untreated coverages, so the same floor admits roughly twice as much there.
```

## Clipping

`--clip-below` and `--clip-above` bound the normalized reactivity. Neither is applied unless given, and neither touches the error, so a signal-to-noise ratio derived from the two shifts with the clip. A missing value is left missing rather than raised to a bound.

## CLI Options

<!-- BEGIN GENERATED cmuts-norm OPTIONS -->
### Arguments

| Argument | Description |
| --- | --- |
| `RATES...` | the reactivities to normalize |

### Input and output

| Option | Description |
| --- | --- |
| `-o, --output HDF5` | write one input's results to this file; repeat once per input (required) |
| `--overwrite` | replace the output files if they already exist |

### Normalization

| Option | Description |
| --- | --- |
| `--norm SCHEME` | how the scale is taken from the rates (ubr, outlier; default ubr) |
| `--min-coverage N` | coverage a position needs before its rate sets the scale (ubr only) (default 500) |

### Clipping

| Option | Description |
| --- | --- |
| `--clip-below N` | raise a normalized reactivity up to this value (default: none) |
| `--clip-above N` | lower a normalized reactivity down to this value (default: none) |

### Information

| Option | Description |
| --- | --- |
| `-h, --help` | show this help and exit |
| `-V, --version` | show the version and exit |

### Advanced

Accepted, and left out of `--help`.

| Option | Description |
| --- | --- |
| `--dump-options` | describe every argument as JSON and exit |
<!-- END GENERATED cmuts-norm OPTIONS -->
