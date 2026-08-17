# cmuts norm

## Purpose

Experiment-independent normalization of reactivity rates.

## Requires

- One or more sets of reactivity rates, in cmuts-compatible HDF5 files

## Usage

Every input given to one run shares a single scale, and each is written to an output of its own. `--output` is repeated once per input and paired with them in order.

```sh
cmuts norm -o apo-normalized.h5 -o holo-normalized.h5 apo.h5 holo.h5
```

Running the inputs separately gives each its own scale instead.

```sh
cmuts norm -o apo-normalized.h5 apo.h5
cmuts norm -o holo-normlized.h5 holo.h5
```

Pool the inputs whenever the rates are to be compared across experiments, since a scale of its own puts each experiment on a different footing.

## Normalization Schemes

`--norm ubr` (the default) takes the 90th percentile of the pooled rates, counting only positions whose coverage exceeds `--min-coverage`. `--norm outlier` drops the highest 2% of the pooled rates as outliers and averages what lies between there and the highest 10%; it reads no coverage, so `--min-coverage` does not apply to it.

A scale that comes out as zero, negative, or undefined is not applied and is recorded as NaN.

## Output

<!-- BEGIN GENERATED cmuts-norm FIELDS -->
{.field}
### `coverage`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `zero`

The number of reads in which this base was present, weighted by PHRED scores.

{.field}
### `reactivity`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

The mutation rate divided by the scale this file records, so a rate reads against the scale rather than as a raw frequency.

{.field}
### `error`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

Standard error of the reactivity values. Purely the statistical error introduced by finite read depths; does not account for experimental or systemic errors.

{.field}
### `reads/lengths`

**Shape** `(n, 2l)` · **Type** `uint64` · **Fill** `zero`

The number of reads passing all filters, binned by length.

{.field}
### `reads/counted`

**Shape** `(n,)` · **Type** `uint64` · **Fill** `zero`

The number of reads passing all filters.

{.field}
### `reads/rejected`

**Shape** `(n,)` · **Type** `uint64` · **Fill** `zero`

The number of reads rejected by at least one filter, or which couldn't be modelled by the HMM

{.field}
### `reads/unmapped`

**Shape** `()` · **Type** `uint64` · **Fill** `zero`

The number of reads not aligned to any reference.

{.field}
### `norm`

**Shape** `()` · **Type** `float32` · **Fill** `NaN`

The scale every rate in this file was divided by.
<!-- END GENERATED cmuts-norm FIELDS -->

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
| `--dump-layout` | describe the output format as JSON and exit |
<!-- END GENERATED cmuts-norm OPTIONS -->
