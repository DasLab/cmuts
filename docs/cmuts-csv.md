# `cmuts csv`

## Purpose

Convert an output into a table that a spreadsheet or a dataframe library reads.

## Requires

- An [output](format.md)

## Usage

Pass the output. The table is written to standard output, one reference after another, in the order the rows of the file follow.

```sh
cmuts csv reactivity.h5 > reactivity.csv
```

The file holds no reference names, so the `reference` column counts the rows from one. Give a FASTA to name them instead.

```sh
cmuts csv -f references.fasta reactivity.h5 > reactivity.csv
```

The FASTA must be the one the input was counted against. It is refused unless it holds one record for each row, of the same length as that row.

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
| [`coverage`](format.md#coverage) | required |
| [`sequence`](format.md#sequence) | required |

All other datasets in an input are ignored.
<!-- END GENERATED cmuts-csv DATASETS -->

## Output

Each row holds one position of one reference. The first three columns identify the position, and each column after them is one dataset of the input, under the name it has in the [format](format.md).

| Column | Meaning |
| --- | --- |
| `reference` | the row number, or the name of the record in the FASTA |
| `position` | the position in the reference, counting from one |
| `base` | the base at that position, decoded from the sequence |

A dataset the input does not hold gets no column. A value the run did not measure is not a number, and is written as an empty field.

The [`sequence`](format.md#sequence) dataset supplies the `base` column. Every row is as wide as the longest reference in the file, and the columns past the end of a reference hold a token of their own. The sequence therefore also gives the length of each reference, and a column past the end gets no row.

## CLI Options

<!-- BEGIN GENERATED cmuts-csv OPTIONS -->
### Arguments

| Argument | Description |
| --- | --- |
| `HDF5` | the cmuts output to convert |

### Input

| Option | Description |
| --- | --- |
| `-f, --fasta FASTA` | name the references from this file, in the order of the rows (default: the row number) |

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
