# cmuts score

## Purpose

Measure a set of reactivities against a structure that is already known.

## Requires

- A cmuts-compatible HDF5 file of reactivities
- The FASTA the file was counted against
- The pairing of one or more of its references, as dot bracket records

## Usage

An output file holds no reference names, and its rows follow the FASTA that produced it, so the same FASTA names them here. Each reference a structure is held for is written as one row of a table.

```sh
cmuts score -f ref.fasta -s structures.db reactivity.h5
```

The structures are read from a file of dot bracket records, matched to the references by name. A record may carry its sequence, which is checked against the reference.

```
>rool120
((((....))))...
>another
GGCAUUAAGCCU
((((....))))
```

A bracket marks a paired base and a dot an unpaired one. Any other character, such as a dash, marks a base the structure does not resolve, and it is left out of the scoring.

One structure can be given in place of a file, which needs the reference it belongs to.

```sh
cmuts score -f ref.fasta -r rool120 --structure '((((....))))...' reactivity.h5
```

## What is Scored

A reagent reports only on the bases it reacts with, so `--bases` restricts the scoring to them. DMS reads adenine and cytosine, and a SHAPE reagent reads all four.

```sh
cmuts score -f ref.fasta -s structures.db -b AC dms.h5
```

`--min-coverage` leaves out the positions too few reads cover to say anything about.

Score a treated sample against its untreated control, not a raw count: `cmuts sub` takes the background off first, and what remains is the signal the treatment added.

```sh
cmuts sub -o subtracted.h5 treated.h5 untreated.h5
cmuts score -f ref.fasta -s structures.db subtracted.h5
```

## Output

One row per reference, as a tab separated table on standard output. A reference whose scored positions are all paired or all unpaired is left out, since a ranking needs both.

| Column | Meaning |
| --- | --- |
| `reference` | the name the FASTA gives it |
| `paired` | scored positions the structure pairs |
| `unpaired` | scored positions the structure leaves open |
| `auroc` | the chance a reactivity ranks an unpaired base above a paired one |
| `auprc` | average precision, taking the unpaired bases as what is sought |
| `mean_paired` | mean reactivity of the paired positions |
| `mean_unpaired` | mean reactivity of the unpaired positions |

Each row covers one reference alone. Reactivity carries the depth and the scale of the sample it came from, so values from two references do not compare, and neither do their rows. Average the rows or resample over them downstream, as the question needs.

## CLI Options

<!-- BEGIN GENERATED cmuts-score OPTIONS -->
### Arguments

| Argument | Description |
| --- | --- |
| `HDF5` | the output being scored |

### Input

| Option | Description |
| --- | --- |
| `-f, --fasta FASTA` | the references the input was counted against, which name its rows (required) |
| `-s, --structures FILE` | dot bracket records, matched to the references by name |
| `--structure DOTBRACKET` | one structure given in place of a file; needs --reference |

### Scoring

| Option | Description |
| --- | --- |
| `-r, --reference NAME` | score this reference alone (default: every reference a structure is held for) |
| `-b, --bases BASES` | score only these bases, as the reagent reports on them (default: every base) |
| `--min-coverage D` | reads a position needs before it is scored (default 0) |

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
<!-- END GENERATED cmuts-score OPTIONS -->
