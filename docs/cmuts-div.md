# `cmuts div`

## Purpose

Normalization of reactivity rates against a denatured control.

## Requires

- Reactivity rates, usually the output of [`cmuts sub`](cmuts-sub.md)
- Denatured control reactivity rates

Both must be in [`cmuts`-compatible HDF5 files](format.md).

## Usage

The control is applied after the background is subtracted.

```sh
cmuts sub -o difference.h5 treated.h5 untreated.h5
cmuts div -o normalized.h5 difference.h5 denatured.h5
```

Each rate is divided by the control's rate at the same position. The result is NaN if either input is NaN there, or if the control's rate is not above zero.

If any dataset does not match the expected shape or the sequences do not match, the program exits early.

## Output

<!-- BEGIN GENERATED cmuts-div DATASETS -->
| Dataset | Source |
| --- | --- |
| [`mismatches/rate`](format.md#mismatchesrate) | Sample over control. |
| [`mismatches/error`](format.md#mismatcheserror) | Propagated from the inputs. |
| [`insertions/rate`](format.md#insertionsrate) | Sample over control. |
| [`insertions/error`](format.md#insertionserror) | Propagated from the inputs. |
| [`deletions/rate`](format.md#deletionsrate) | Sample over control. |
| [`deletions/error`](format.md#deletionserror) | Propagated from the inputs. |
| [`coverage`](format.md#coverage) | Summed over the inputs. |
| [`sequence`](format.md#sequence) | Copied from the inputs. |
| [`reads/lengths`](format.md#readslengths) | Summed over the inputs. |
| [`reads/counted`](format.md#readscounted) | Summed over the inputs. |
| [`reads/rejected`](format.md#readsrejected) | Summed over the inputs. |
| [`reads/unmapped`](format.md#readsunmapped) | Summed over the inputs. |
<!-- END GENERATED cmuts-div DATASETS -->

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
| `--dump-layout` | describe the output format as JSON and exit |
<!-- END GENERATED cmuts-div OPTIONS -->
