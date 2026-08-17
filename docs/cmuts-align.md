# cmuts align

## Purpose

Aligning raw sequencing data against the reference library.

## Requires

- One or more FASTQ files containing the sequenced reads
- The FASTA library

As well as both `minimap2` and `samtools` on the path, and `fastp` for paired-end input.

## Usage

One invocation handles one sample. Pass a single FASTQ file for single-end reads,

```sh
cmuts align -f references.fasta -o treated.bam -x map-ont treated.fastq
```

or two for paired-end reads.

```sh
cmuts align -f references.fasta -o treated.bam -x sr treated_R1.fastq.gz treated_R2.fastq.gz
```

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

## Options

### Arguments

| Argument | Description |
| --- | --- |
| `READS` | reads to align |
| `MATE` | the second file of a pair, merged with the first |

### Input and output

| Option | Description |
| --- | --- |
| `-f, --fasta FASTA` | reference sequences, or a minimap2 index (required) |
| `-o, --output BAM` | write sorted alignments to this file (required) |
| `--overwrite` | replace the output file if it already exists |

### Alignment

| Option | Description |
| --- | --- |
| `-x, --preset PRESET` | how minimap2 aligns (sr\|map-ont\|map-hifi\|map-pb) (required) |

### Performance

| Option | Description |
| --- | --- |
| `-t, --threads N` | threads for merging, alignment and sorting (default 1) |

### Information

| Option | Description |
| --- | --- |
| `-h, --help` | show this help and exit |
| `-V, --version` | show the version and exit |
