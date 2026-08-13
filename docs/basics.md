# Using cmuts

Running cmuts requires one or more coordinate-sorted alignment files (SAM, BAM, or CRAM) and the reference FASTA to which they were aligned. For a full list of the arguments each command takes, please read their respective pages.

## Standard Usage

The basic way to use cmuts on a standard MaP-seq dataset involves passing the treated and untreated samples through the HMM in order to get per-experiment reactivities,

```sh
cmuts-hmm -f references.fasta -o treated.h5 treated.bam
cmuts-hmm -f references.fasta -o untreated.h5 untreated.bam
```

and then performing background subtraction.

```sh
cmuts-sub -o combined.h5 treated.h5 untreated.h5
```

All HDF5 files have the same format, where `n` is the number of references and `l` the length of the longest of them.

<!-- BEGIN GENERATED cmuts-hmm LAYOUT -->
| Dataset | Shape | Type | Fill |
| --- | --- | --- | --- |
| `coverage` | `(n, l)` | `float32` | `zero` |
| `reactivity` | `(n, l)` | `float32` | `NaN` |
| `error` | `(n, l)` | `float32` | `NaN` |
| `reads/lengths` | `(n, 2l)` | `uint64` | `zero` |
| `reads/counted` | `(n,)` | `uint64` | `zero` |
| `reads/rejected` | `(n,)` | `uint64` | `zero` |
| `reads/unmapped` | `()` | `uint64` | `zero` |
<!-- END GENERATED cmuts-hmm LAYOUT -->

See the [output](output.md) page for more detail on what each dataset contains.

## Split Alignments

If reads are split across lanes, then multiple alignment files can be passed to the HMM, where they are treated as if they were one large alignment file.

```sh
cmuts-hmm -f references.fasta -o counts.h5 lane1.bam lane2.bam lane3.bam
```

## Denatured Controls

To account for a denatured control, first use the HMM to get its mutation rates alongside the treated and untreated experiments, and then pass it to the `-d` flag of `cmuts-sub`.

```sh
cmuts-hmm -f references.fasta -o treated.h5 treated.bam
cmuts-hmm -f references.fasta -o untreated.h5 untreated.bam
cmuts-hmm -f references.fasta -o denatured.h5 denatured.bam
cmuts-sub -o combined.h5 -d denatured.h5 treated.h5 untreated.h5
```
