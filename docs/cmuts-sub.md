# cmuts-sub

## Purpose

Background subtraction of reactivity rates.

## Requires

- Treated reactivity rates
- Untreated reactivity rates

Both must be in cmuts-compatible HDF5 files.

## Error and Coverage

The error of the background-subtracted rates is computed using the standard quadrature formula. The coverage and the four read-related datasets simply take the sum of the treated and untreated values.

To normalize the result against a denatured control, pass it to [cmuts-div](cmuts-div.md).

```{warning}
High coverage in the output dataset does not imply high-quality data, since it is insensitive to imbalances in the treated and untreated experiments. Either ensure each experiment separately has high coverage or use the signal-to-noise ratio as a more robust quality metric.
```

## Output

<!-- BEGIN GENERATED cmuts-sub FIELDS -->
{.field}
### `coverage`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `zero`

The number of reads in which this base was present, weighted by PHRED scores.

{.field}
### `reactivity`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

The mutation rate of the treated sample less that of the untreated one, so what remains is the signal the treatment added.

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
<!-- END GENERATED cmuts-sub FIELDS -->

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
| `--clip` | raise a negative reactivity to zero |

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
<!-- END GENERATED cmuts-sub OPTIONS -->
