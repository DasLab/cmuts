# `cmuts norm`

## Purpose

Experiment-independent normalization of reactivity rates.

## Requires

- One or more sets of reactivity rates, in [`cmuts`-compatible HDF5 files](format.md)

## Usage

Every input given to one run shares a single norm, and each is written to an output of its own. `-o` is repeated once per input and paired with them in order.

```sh
cmuts norm -o apo-normalized.h5 -o holo-normalized.h5 apo.h5 holo.h5
```

Running the inputs separately gives each its own norm instead.

```sh
cmuts norm -o apo-normalized.h5 apo.h5
cmuts norm -o holo-normlized.h5 holo.h5
```

This is only suggested if the two experiments had significantly different conditions which are not directly comparable.

## Normalization Schemes

`--norm ubr` (the default) takes the 90th percentile of the pooled rates, counting only positions whose coverage exceeds `--min-coverage`. `--norm outlier` drops the highest 2% of the pooled rates as outliers and averages what lies between there and the highest 10%; it reads no coverage, so `--min-coverage` does not apply to it.

A norm that comes out as zero, negative, or undefined is not applied and is recorded as NaN.

## Output

<!-- BEGIN GENERATED cmuts-norm DATASETS -->
| Dataset | Source |
| --- | --- |
| [`mismatches/rate`](format.md#mismatchesrate) | Divided by the norm. |
| [`mismatches/error`](format.md#mismatcheserror) | Divided by the norm. |
| [`norm`](format.md#norm) | Estimated per the specified scheme. |
| [`coverage`](format.md#coverage) | Copied from the input. |
| [`sequence`](format.md#sequence) | Copied from the input. |
| [`reads/lengths`](format.md#readslengths) | Copied from the input. |
| [`reads/counted`](format.md#readscounted) | Copied from the input. |
| [`reads/rejected`](format.md#readsrejected) | Copied from the input. |
| [`reads/unmapped`](format.md#readsunmapped) | Copied from the input. |
<!-- END GENERATED cmuts-norm DATASETS -->

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
| `--norm SCHEME` | how the norm is taken from the rates (ubr, outlier; default ubr) |
| `--min-coverage N` | coverage a position needs before its rate sets the norm (ubr only) (default 500) |

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
