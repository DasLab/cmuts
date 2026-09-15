# `cmuts align`

## Purpose

Aligning raw sequencing data against a reference library.

## Requires

- The sequenced reads, either in FASTQ or unaligned BAM format
- The FASTA library

As well as both `minimap2` and `samtools` on the path, and `fastp` if processing paired-end input.

## Usage

One invocation handles one sample. Pass a single FASTQ file for single-end reads,

```sh
cmuts align -f references.fasta -o treated.bam -x map-ont treated.fastq
```

or two for paired-end reads.

```sh
cmuts align -f references.fasta -o treated.bam -x sr treated_R1.fastq.gz treated_R2.fastq.gz
```

PacBio and some nanopore instruments deliver reads as an unaligned BAM rather than a FASTQ, which is passed the same way.

```sh
cmuts align -f references.fasta -o treated.bam -x map-hifi treated.hifi_reads.bam
```

## Presets

The required `-x` argument names the platform the reads came from, which is passed on to `minimap2` for alignment purposes. It must be one of the following.

| Preset | Platform |
| --- | --- |
| `sr` | Illumina, MGI, Complete Genomics |
| `map-ont` | Oxford Nanopore |
| `map-hifi` | PacBio HiFi |
| `map-pb` | PacBio CLR |
| `map-iclr` | Illumina Complete Long Reads |
| `lr:hq` | Any long read platform with an error rate below one percent |

Paired-end inputs are refused unless the preset is `sr`.

## Unaligned BAM

The program reads only the name, the bases, and the base qualities from an unaligned BAM. Per-base tags such as the PacBio `ip` and `pw` fields, or the base modification calls in `MM` and `ML`, are lost. If a record has no sequence, has no base qualities, or is marked as secondary, supplementary, or reverse, the program exits and no output is written.

Paired-end data in unaligned BAM format is not supported. BAM files which are already aligned are refused.

## Checksums

The output header contains the MD5 checksum of each reference in the `M5` field of its `@SQ` line. The checksum is used by [`cmuts hmm`](cmuts-hmm.md) to ensure alignment and processing use the same FASTA.

## Merging

Paired-end input is merged with `fastp` before alignment. A pair whose mates do not overlap cannot be merged and is discarded. `fastp` reports how many pairs it merged on standard error.

## Sorting

The alignments are sorted by `samtools sort`, with `-t` and `-m` controlling the number of threads and the maximum memory per thread.

## CLI Options

<!-- BEGIN GENERATED cmuts-align OPTIONS -->
### Arguments

| Argument | Description |
| --- | --- |
| `READS` | reads to align, as FASTQ or unaligned BAM |
| `MATE` | the second file of a pair, merged with the first |

### Input and output

| Option | Description |
| --- | --- |
| `-f, --fasta FASTA` | reference sequences (required) |
| `-o, --output BAM` | write sorted alignments to this file (required) |
| `--overwrite` | replace the output file if it already exists |

### Alignment

| Option | Description |
| --- | --- |
| `-x, --preset PRESET` | minimap2 preset for the sequencing technology (sr, map-ont, map-hifi, map-pb, map-iclr, lr:hq; required) |

### Performance

| Option | Description |
| --- | --- |
| `-t, --threads N` | threads for merging, alignment and sorting (default 1) |
| `-m, --memory MiB` | memory each sorting thread holds before it writes to disk (default 768) |

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
<!-- END GENERATED cmuts-align OPTIONS -->
