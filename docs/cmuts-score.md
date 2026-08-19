# cmuts score

## Purpose

Measure reactivities against a known structure.

## Requires

- A cmuts-compatible HDF5 file of reactivities
- The FASTA that the file was counted against
- The pairing of one or more of its references, as dot bracket records

## Usage

An output file stores no reference names. Its rows are in the order of the FASTA that produced it, so `cmuts score` takes the names and the sequences from that FASTA. Each reference with a structure gives one row of the output table.

```sh
cmuts score -f ref.fasta -s structures.db reactivity.h5
```

The structures file holds dot bracket records. Each record is matched to a reference by name. A record can also hold the sequence.

```
>rool120
((((....))))...
>another
GGCAUUAAGCCU
((((....))))
```

A bracket marks a paired base, and a dot marks an unpaired base. Any other character, such as a dash, marks a base that the structure does not resolve. Those positions are not scored.

If a record holds a sequence, it is compared with the reference, and the run stops at the first base that differs. U and T are equivalent in this comparison, so an RNA structure matches a DNA reference.

The file does not have to cover the whole library. A reference without a record is skipped, and a record whose name is not in the library is ignored. The rows are in the order of the FASTA, whatever order the records are in.

## What is Scored

Each reagent modifies only some bases, so `--bases` limits the scoring to those bases. DMS modifies adenine and cytosine. A SHAPE reagent modifies all four bases.

```sh
cmuts score -f ref.fasta -s structures.db -b AC dms.h5
```

`--min-coverage` removes the positions that too few reads cover.

Score a treated sample against its untreated control, and not a raw count. `cmuts sub` removes the background first, and the result is the signal that the treatment added.

```sh
cmuts sub -o subtracted.h5 treated.h5 untreated.h5
cmuts score -f ref.fasta -s structures.db subtracted.h5
```

## Output

One row for each reference, as a tab separated table on standard output. A reference is not scored if all of its scored positions are paired, or if all of them are unpaired, because a ranking needs both classes.

| Column | Meaning |
| --- | --- |
| `reference` | the name in the FASTA |
| `paired` | scored positions that the structure pairs |
| `unpaired` | scored positions that the structure leaves unpaired |
| `auroc` | the chance that an unpaired base has a higher reactivity than a paired base |
| `auprc` | average precision, with the unpaired bases as the positive class |
| `mean_paired` | mean reactivity of the paired positions |
| `mean_unpaired` | mean reactivity of the unpaired positions |

Each row covers one reference. The reactivity of a position depends on the read depth and the scale of its sample, so values from two references are not comparable. Average the rows, or resample over them, in a later step.

## CLI Options

<!-- BEGIN GENERATED cmuts-score OPTIONS -->
### Arguments

| Argument | Description |
| --- | --- |
| `HDF5` | the output being scored |

### Input

| Option | Description |
| --- | --- |
| `-f, --fasta FASTA` | the references the input was counted against, in the order of its rows (required) |
| `-s, --structures FILE` | dot bracket records, matched to the references by name (required) |

### Scoring

| Option | Description |
| --- | --- |
| `-b, --bases BASES` | score only the bases the reagent modifies (default: every base) |
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
