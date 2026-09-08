# `cmuts align`

## Purpose

Aligning raw sequencing data against the reference library.

## Requires

- The sequenced reads, as one or more FASTQ files or as one unaligned BAM
- The FASTA library

As well as both `minimap2` and `samtools` on the path, `fastp` for paired-end input, and `gzip` for compressed input.

## Usage

One invocation handles one sample. Pass a single FASTQ file for single-end reads,

```sh
cmuts align -f references.fasta -o treated.bam -x map-ont treated.fastq
```

or two for paired-end reads.

```sh
cmuts align -f references.fasta -o treated.bam -x sr treated_R1.fastq.gz treated_R2.fastq.gz
```

## Unaligned BAM

PacBio and nanopore instruments deliver reads as an unaligned BAM rather than a FASTQ, and one may be passed in place of one.

```sh
cmuts align -f references.fasta -o treated.bam -x map-hifi treated.hifi_reads.bam
```

The records are read back as FASTQ through `samtools fastq` and aligned as any other reads are. Nothing beyond the name, the bases, and the base qualities carries into the output, so the per-base tags an instrument writes, such as the PacBio `ip` and `pw` fields or the base modification calls in `MM` and `ML`, are not preserved.

A BAM holds the reads of both mates, so it is passed on its own and never as one file of a pair.

An already-aligned BAM is refused. Alignments are what `cmuts hmm` counts, so a BAM whose header names the references its reads were placed against is a confusion between the two subcommands rather than work for this one.

## Checksums

The output header declares the MD5 checksum of each reference in the `M5` field of its `@SQ` line. `cmuts hmm --verify checksum` compares these against the FASTA it is given.

## Merging

Paired-end input is merged with `fastp` before alignment. A pair whose mates do not overlap cannot be merged and is discarded. `fastp` reports how many pairs it merged on standard error.

## Presets

`-x` is required, and names the platform the reads came from, which is passed on to `minimap2` for alignment purposes. The following lists common choices; consult the `minimap2` documentation for all choices.

| Preset | Platform |
| --- | --- |
| `sr` | Illumina, MGI, Complete Genomics |
| `map-ont` | Oxford Nanopore |
| `map-hifi` | PacBio HiFi |
| `map-pb` | PacBio CLR |

Paired-end inputs are refused unless the preset is `sr`.

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
| `-x, --preset PRESET` | how minimap2 aligns (sr, map-ont, map-hifi, map-pb; required) |

### Performance

| Option | Description |
| --- | --- |
| `-t, --threads N` | threads for merging, alignment and sorting (default 1) |

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
