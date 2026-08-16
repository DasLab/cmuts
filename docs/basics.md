# Basics

The cmuts pipeline comprises five independent programs:

{.field}
[`cmuts-align`](cmuts-align.md)\
**Purpose:** Aligning raw sequencing data against the reference library\
**Requires:** One or more FASTQ files, the FASTA library

{.field}
[`cmuts-hmm`](cmuts-hmm.md)\
**Purpose:** Computing reactivity rates from alignment files via the pair-HMM\
**Requires:** One or more coordinate-sorted alignment (SAM, BAM, or CRAM) files, the FASTA library

{.field}
[`cmuts-sub`](cmuts-sub.md)\
**Purpose:** Background subtraction of reactivity rates\
**Requires:** Treated and untreated reactivity rates, in cmuts-compatible HDF5 files

{.field}
[`cmuts-div`](cmuts-div.md)\
**Purpose:** Normalization of reactivity rates against a denatured control\
**Requires:** Reactivity rates and denatured control rates, in cmuts-compatible HDF5 files

{.field}
[`cmuts-norm`](cmuts-norm.md)\
**Purpose:** Normalization of reactivity rates against a scale taken from the rates themselves\
**Requires:** One or more sets of reactivity rates, in cmuts-compatible HDF5 files

This page goes over basic, end-to-end usage of these programs on standard data. For a full list of the arguments each command takes, please read their respective pages.

## Standard Usage

The canonical use case of the cmuts pipeline is to generate reactivity profiles from a MaP-seq experiment, where the cDNA of treated and untreated RNA has been sequenced.

The first step is to align the reads against the reference library.

```sh
cmuts-align -f references.fasta -x sr -o treated.bam treated.fastq.gz
cmuts-align -f references.fasta -x sr -o untreated.bam untreated.fastq.gz
```

`-x` names the platform the reads were sequenced with, and is required. See the [cmuts-align](cmuts-align.md) page for the presets it accepts and the rest of its options.

Then, pass the alignments to the HMM in order to compute reactivity rates.

```sh
cmuts-hmm -f references.fasta -o treated.h5 treated.bam
cmuts-hmm -f references.fasta -o untreated.h5 untreated.bam
```

With this done, perform background subtraction.

```sh
cmuts-sub -o reactivity.h5 treated.h5 untreated.h5
```

The final step is normalizing the reactivity.

```sh
cmuts-norm -o normalized-reactivity.h5 reactivity.h5
```

All HDF5 files in cmuts have the same format, where `n` is the number of references and `l` the length of the longest of them.

<!-- BEGIN GENERATED cmuts-norm LAYOUT -->
| Dataset | Shape | Type | Fill |
| --- | --- | --- | --- |
| `coverage` | `(n, l)` | `float32` | `zero` |
| `reactivity` | `(n, l)` | `float32` | `NaN` |
| `error` | `(n, l)` | `float32` | `NaN` |
| `reads/lengths` | `(n, 2l)` | `uint64` | `zero` |
| `reads/counted` | `(n,)` | `uint64` | `zero` |
| `reads/rejected` | `(n,)` | `uint64` | `zero` |
| `reads/unmapped` | `()` | `uint64` | `zero` |
| `norm` | `()` | `float32` | `NaN` |
<!-- END GENERATED cmuts-norm LAYOUT -->

See the [output](output.md) page for more detail on what each dataset contains.

## Pre-Aligned Data

Skip running `cmuts-align`. Ensure your data is sorted, which can be done with `samtools sort`.

## Split Alignments

Multiple alignment files can be passed to the HMM, where they are treated as if they were one large alignment file.

```sh
cmuts-hmm -f references.fasta -o counts.h5 lane1.bam lane2.bam lane3.bam
```

## Denatured Control

To account for a denatured control, first use the HMM to get its mutation rates alongside the treated and untreated experiments,

```sh
cmuts-hmm -f references.fasta -o denatured.h5 denatured.bam
```

and then divide the background-subtracted rates by it.

```sh
cmuts-div -o normalized.h5 combined.h5 denatured.h5
```

## No Control

With no control, skip background subtraction and use the output of `cmuts-hmm` directly.

## Correlation-Based Data

To analyze M2-seq, RING-MaP, or MOHCA-seq data, pass the `--pairwise` flag to `cmuts-hmm`, which causes it to compute and store pairwise mutation correlations.
