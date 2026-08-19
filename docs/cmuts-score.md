# cmuts score

## Purpose

Score reactivity rates against known structures.

## Requires

- Reactivity rates, in a cmuts-compatible HDF5 file
- Dot-bracket secondary structures
- The reference library FASTA

## Usage

Pass the reactivity, structures, and reference library to `cmuts score`.

```sh
cmuts score -f ref.fasta -s structures.db reactivity.h5
```

The structure file holds dot-bracket records, each matched to its reactivity profile by name.

```
>sequence1
((((....))))...
>sequence2
((..(((...)))))
```

Positions containing neither dots nor brackets are not scored.

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
| `-b, --bases BASES` | score only at these bases (A, C, G, U; default A,C,G,U) |
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
