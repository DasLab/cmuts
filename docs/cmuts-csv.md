# `cmuts csv`

## Purpose

Converting a `cmuts` HDF5 file into a comma-separated table.

## Requires

- A [`cmuts`-compatible HDF5 file](format.md)

## Usage

Provide the HDF5 file and redirect the output to a file.

```sh
cmuts csv reactivity.h5 > reactivity.csv
```

Providing a FASTA populates the `reference` column with the corresponding name.

```sh
cmuts csv -f references.fasta reactivity.h5 > reactivity.csv
```

## Input

<!-- BEGIN GENERATED cmuts-csv DATASETS -->
| Dataset | Input |
| --- | --- |
| [`mismatches/rate`](format.md#mismatchesrate) | if present |
| [`mismatches/error`](format.md#mismatcheserror) | if present |
| [`insertions/rate`](format.md#insertionsrate) | if present |
| [`insertions/error`](format.md#insertionserror) | if present |
| [`deletions/rate`](format.md#deletionsrate) | if present |
| [`deletions/error`](format.md#deletionserror) | if present |
| [`terminations/rate`](format.md#terminationsrate) | if present |
| [`terminations/error`](format.md#terminationserror) | if present |
| [`coverage`](format.md#coverage) | required |
| [`sequence`](format.md#sequence) | required |

All other datasets in an input are ignored.
<!-- END GENERATED cmuts-csv DATASETS -->

## Output

Each row is one position in one reference. The first three columns identify the position as follows.

| Column | Meaning |
| --- | --- |
| `reference` | the reference number, or the name of the record in the FASTA |
| `position` | the position in the reference, counting from one |
| `base` | the base at that position, taken from the [`sequence`](format.md#sequence) dataset |

Each column after these is one dataset of the input, under the name it has in the [format](format.md). NaN values are converted to empty fields.

In ragged libraries references are not padded.

## CLI Options

<!-- BEGIN GENERATED cmuts-csv OPTIONS -->
### Arguments

| Argument | Description |
| --- | --- |
| `HDF5` | the cmuts output to convert |

### Input

| Option | Description |
| --- | --- |
| `-f, --fasta FASTA` | populate the reference field with the names from this file (default: the reference number) |

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
| `--dump-layout` | describe the input format as JSON and exit |
<!-- END GENERATED cmuts-csv OPTIONS -->
