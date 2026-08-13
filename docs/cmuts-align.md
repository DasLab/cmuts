# cmuts-align

`cmuts-align` aligns reads to a reference with [minimap2](https://github.com/lh3/minimap2) and converts this to a coordinate-sorted BAM suitable for `cmuts-hmm` using [samtools](https://github.com/samtools/samtools). It needs both `minimap2` and `samtools` on the path.

One invocation handles one sample. Pass a single FASTQ file for single-end reads,

```sh
cmuts-align -f references.fasta -o treated.bam -x map-ont treated.fastq
```

or two for paired-end reads.

```sh
cmuts-align -f references.fasta -o treated.bam -x sr treated_R1.fastq.gz treated_R2.fastq.gz
```

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
| `MATE` | the second file of a pair |

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
| `-t, --threads N` | threads for alignment and sorting (default 1) |

### Information

| Option | Description |
| --- | --- |
| `-h, --help` | show this help and exit |
| `-V, --version` | show the version and exit |
