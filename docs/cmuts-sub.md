# cmuts-sub

## Purpose

Background subtraction of reactivity rates.

## Requires

- Treated reactivity rates
- Untreated reactivity rates

Both must be in cmuts-compatible HDF5 files.

## Error and Coverage

The error of the background-subtracted rates is computed using the standard quadrature formula. The coverage and the four read-related datasets simply take the sum of the treated and untreated values.

```{warning}
High coverage in the output dataset does not imply high-quality data, since it is insensitive to imbalances in the treated and untreated experiments. Either ensure each experiment separately has high coverage or use the signal-to-noise ratio as a more robust quality metric.
```

## CLI Options

<!-- BEGIN GENERATED cmuts-sub OPTIONS -->
### Arguments

| Argument | Description |
| --- | --- |
| `TREATED` | the modified sample |
| `UNTREATED` | the background |

### Input and output

| Option | Description |
| --- | --- |
| `-o, --output HDF5` | write results to this file (required) |
| `--overwrite` | replace the output file if it already exists |

### Subtraction

| Option | Description |
| --- | --- |
| `-d, --denatured HDF5` | normalize against a denatured control (default: none) |
| `--clip` | raise a negative reactivity to zero |

### Information

| Option | Description |
| --- | --- |
| `-h, --help` | show this help and exit |
| `-V, --version` | show the version and exit |
<!-- END GENERATED cmuts-sub OPTIONS -->
