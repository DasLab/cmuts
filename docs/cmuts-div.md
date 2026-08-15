# cmuts-div

## Purpose

Normalization of reactivity rates against a denatured control.

## Requires

- Reactivity rates, usually the output of [cmuts-sub](cmuts-sub.md)
- Denatured control reactivity rates

Both must be in cmuts-compatible HDF5 files.

## Usage

The control is applied after the background is subtracted.

```sh
cmuts-sub -o difference.h5 treated.h5 untreated.h5
cmuts-div -o normalized.h5 difference.h5 denatured.h5
```

Each rate is divided by the control's rate at the same position. The result is NaN if either input is NaN there, or if the control's rate is not above zero.

## Error and Coverage

The error is computed via the standard propagation of error for a ratio.  The coverage and the four read-related datasets simply take the sum of the treated and untreated values.

```{warning}
High coverage in the output dataset does not imply high-quality data, since it is insensitive to imbalances between the experiments. Either ensure each experiment separately has high coverage or use the signal-to-noise ratio as a more robust quality metric.
```

## CLI Options

<!-- BEGIN GENERATED cmuts-div OPTIONS -->
### Arguments

| Argument | Description |
| --- | --- |
| `RATES` | the reactivities to normalize |
| `CONTROL` | the denatured control |

### Input and output

| Option | Description |
| --- | --- |
| `-o, --output HDF5` | write results to this file (required) |
| `--overwrite` | replace the output file if it already exists |

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
<!-- END GENERATED cmuts-div OPTIONS -->
