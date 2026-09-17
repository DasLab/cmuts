# Changelog

## [2.0.0] - 2026-09-16

cmuts v2 is a rewrite of cmuts. The pipeline is one binary. The counting runs on a pair HMM, and each output file holds finished rates rather than raw counts.

**Pair HMM**: A pair HMM marginalizes over the alignments in a band around the CIGAR, in place of v1's deletion spreading. It handles ambiguous insertions, which v1 did not. `--band` sets how far either side of the CIGAR it looks.

**One binary**: The pipeline is now subcommands of one `cmuts` binary: `cmuts align`, `cmuts hmm`, `cmuts sub`, `cmuts div`, `cmuts norm`, `cmuts csv`, `cmuts plot`, and `cmuts gen`. Each release carries this binary statically linked for Linux (x86_64, aarch64) and macOS (arm64). Running it needs no compiler and no libraries.

**Better parallelism**: Threads replace `MPI` processes. Several threads can count the reads on one reference, where v1 gave each reference to one process. A reference with no reads is skipped.

**Fewer dependencies**: `OMP`, `MPI`, `HDF5-MPI`, and `htscodecs` are all no longer needed. There is no separate parallel build to manage. The build system is reduced from `cmake` to `make`.

**Minimal Python**: The Python package is gone. The pipeline is one C binary, alongside one Python script that serves the interactive report. Python is otherwise only in the test suite.

### Commands

| v1 | v2 |
| --- | --- |
| `cmuts align` | `cmuts align` |
| `cmuts core` | `cmuts hmm` |
| `cmuts normalize` | `cmuts sub` (background subtraction), `cmuts div` (denatured control), and `cmuts norm` (normalization) |
| `cmuts generate` | `cmuts gen` |
| `cmuts plot` | `cmuts plot` |
| `cmuts visualize` | [cif-overlay](https://github.com/hmblair/cif-overlay), a program of its own |
| `cmuts test` | `make check`, from a checkout |

### Added

- `cmuts div` divides reactivity rates by a denatured control.
- `cmuts csv` writes an output as a table of comma separated values, one row per position of each reference.
- `cmuts align` may read an unaligned BAM in place of a FASTQ.
- `--pairwise` names the statistics to write, either `correlation` or `conditional`. v1 wrote raw joint counts for `cmuts normalize` to process.
- `--params` reads the pair-HMM rates from a file, and `--dump-params` writes the defaults in the same form.
- `--alignment-type` allows filtering primary or supplementary alignments, with the latter rejected by default.
- `--verify` checks the FASTA against the alignment header.
- `--min-coverage` sets the coverage a position needs before its rate sets the normalization factor.
- `--norm value` divides by the number `--value` gives, in place of computing one from the rates. It applies the normalization factor of one run to another.
- The `reads` datasets hold the number of reads counted and rejected at each reference, under `reads/primary`, and the same two counts for the further pieces of split reads, under `reads/supplementary`. They also hold the lengths of the counted reads, and the number of reads in the alignment that map to no reference.

### Changed

- Each kind of difference has its own rate and error: `mismatches`, `insertions`, `deletions`, and `terminations`. v1 wrote one `reactivity` dataset.
- The insertion and deletion rates divide by the reads that cover the position and continue past it. v1 divided by the coverage.
- Each base's PHRED score weights its contribution. `--min-phred` marks a base below it as carrying no information, where v1 rejected the base.
- `--min-mapq` defaults to 20, where v1 defaulted to 10.
- `--nan-5p` and `--nan-3p` are v1's `--blank-5p` and `--blank-3p`, applied by `cmuts hmm` rather than at normalization.
- `--min-depth` is v1's `--blank-cutoff`, and defaults to 1 rather than 10.
- `--strand` names the strands to keep, in place of `--no-reverse` and `--only-reverse`.
- `-j` sets how many threads count reads, in place of `--threads`, which ran that many `MPI` processes.
- `cmuts sub` clips a negative reactivity to zero unless `--keep-negative` is given, in place of `--clip-below` and `--clip-above`.
- `cmuts align` uses `minimap2` instead of `bowtie2`, gaining presets `sr`, `map-ont`, `map-hifi`, `map-pb`, `map-iclr`, and `lr:hq`. It merges paired-end mates through `vsearch` before alignment.
- `cmuts hmm` requires coordinate-sorted input and refuses paired reads. Merge the mates before alignment, as `cmuts align` does for paired-end input.
- Several alignment files given to one run are read as one merged alignment, where v1 wrote one group per input file. Replicates merge the same way.
- `cmuts sub`, `cmuts div`, and `cmuts norm` read and write whole HDF5 files, where v1's `--experiment` named datasets inside one counts file.
- Each output file holds one flat layout, in place of v1's groups for each input file, each experiment, and the metadata. The `program` and `version` attributes record what wrote the file.
- `sequence` always holds the reference bases, where v1 wrote them only under `--tokenize`.
- Every dataset is compressed with deflate at level 3, in place of `--compression`.
- `cmuts hmm` streams its input and writes no `.cmix` or `.cmfa` index files beside it.

### Removed

- The `raw`, `sm-dms`, and `sm-shape` normalization schemes, and `--per-experiment-norm` and `--per-reference-norm`. `cmuts norm` offers `ubr` and `outlier`, and takes one factor over every input given to a run.
- The deletion-spreading modes `--uniform-spread`, `--no-spread`, and `--disable-ambiguous`, and `--collapse`.
- The counting toggles `--no-mismatches`, `--no-insertions`, and `--no-deletions`.
- The PHRED filters `--quality-window`, `--no-match-filter`, `--no-insertion-filter`, and `--no-deletion-filter`.
- `--max-hamming`, `--secondary`, `--downsample`, `--ignore-bases`, and `--max-indel-length`.
- The datasets `SNR`, `pairwise-snr`, `mutual-information`, `covariance`, `heatmap`, and `roi-mask`.
- The reference names, which v1 kept under `meta`.
- Demultiplexing and trimming in `cmuts align`, and the ultraplex, cutadapt, and bowtie2 dependencies.
