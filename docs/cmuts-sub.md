# `cmuts sub`

## Purpose

Background subtraction of reactivity rates.

## Requires

- Treated reactivity rates
- Untreated reactivity rates

Both must be in [`cmuts`-compatible HDF5 files](format.md).

## Usage

```sh
cmuts sub -o reactivity.h5 treated.h5 untreated.h5
```

The untreated reactivity is subtracted from the treated reactivity. The result is NaN if either input is NaN there.

If any dataset does not match the expected shape or the sequences do not match, the program exits early.

## Output

<!-- BEGIN GENERATED cmuts-sub DATASETS -->
| Dataset | Source |
| --- | --- |
| [`mismatches/rate`](format.md#mismatchesrate) | Treated less untreated. |
| [`mismatches/error`](format.md#mismatcheserror) | Propagated from the inputs. |
| [`coverage`](format.md#coverage) | Summed over the inputs. |
| [`sequence`](format.md#sequence) | Copied from the inputs. |
| [`reads/lengths`](format.md#readslengths) | Summed over the inputs. |
| [`reads/counted`](format.md#readscounted) | Summed over the inputs. |
| [`reads/rejected`](format.md#readsrejected) | Summed over the inputs. |
| [`reads/unmapped`](format.md#readsunmapped) | Summed over the inputs. |
<!-- END GENERATED cmuts-sub DATASETS -->

```{warning}
High coverage in the output does not imply high-quality data, since it is insensitive to imbalances in the treated and untreated experiments. Either ensure each experiment separately has high coverage or use the signal-to-noise ratio as a more robust quality metric.
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
