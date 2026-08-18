# Changelog

## [2.0.0] - TBD

cmuts v2 is a rewrite of cmuts that makes it more accurate, portable, and, in specific cases, faster. The headline features are

**Pair-HMM**: Gone is the ad-hoc and exponentially-expensive deletion spreading from v1, replaced by a variable-bandwidth pair-HMM that marginalizes over all alignments in a band around the CIGAR. The change also adds handling for ambiguous insertions and merging of nearby mutations for free.

**One binary**: The pipeline is now subcommands of one `cmuts` binary: `cmuts align`, `cmuts hmm`, `cmuts sub`, `cmuts div`, `cmuts norm`, and `cmuts gen`. Each release carries this binary statically linked for Linux (x86_64, aarch64) and macOS (arm64), so running it needs no compiler and no libraries.

**Better parallelism**: Multi-threading is re-implemented via `pthreads`, with an improved model that allows for intra-reference parallelism, speeding up runs with few references but many reads. References with no reads are skipped, where before each occupied an `MPI` process, speeding up runs with many references but few reads.

**Fewer dependencies**: Switching to `pthreads` means `OMP`, `MPI`, and `HDF5-MPI` are all no longer needed, and that there is no separate parallel build to manage. The build system is reduced from `cmake` to `make`, and `htscodecs` is dropped, removing `autoconf`, `automake`, and `libtool` as transitive dependencies too.

**No Python**: The Python package is gone. The pipeline is one C binary and one bash script; Python remains only in the test suite.

### Commands

| v1 | v2 |
| --- | --- |
| `cmuts align` | `cmuts align` |
| `cmuts core` | `cmuts hmm` |
| `cmuts normalize` | `cmuts sub` (background subtraction) and `cmuts norm` (normalization) |
| `cmuts generate` | `cmuts gen` |
| `cmuts plot`, `cmuts visualize` | to be ported before the 2.0.0 release |
| `cmuts test` | `make check`, from a checkout |

### Added

- `cmuts div` divides reactivity rates by a denatured control.
- `cmuts align` aligns through minimap2, so the long-read presets `map-ont`, `map-hifi`, and `map-pb` join short reads, and it merges paired-end mates through fastp before alignment.
- `--pairwise` names the statistics to write. In v1 it was a flag on `cmuts core` that wrote raw joint counts for `cmuts normalize` to process; `cmuts hmm` now writes the finished statistics directly. Mutual information returns as a statistic before the 2.0.0 release.
- `--params` reads the pair-HMM rates from a file, and `--dump-params` writes the defaults in the same form.
- `--verify` checks the FASTA against the alignment header, by each reference's name, length, and MD5 checksum where present.

### Changed

- PHRED scores weight each base's contribution to the counts, in place of the `--min-phred` threshold, `--quality-window`, and the per-type filter toggles.
- `--substitution-weight`, `--deletion-weight`, and `--insertion-weight`, each 0 to 1, set what each kind of difference counts towards the mutation total, in place of the binary `--exclude-mismatches`, `--exclude-deletions`, and `--include-insertions`.
- `--strand` names the strands to keep, in place of `--no-reverse` and `--only-reverse`.
- `cmuts hmm` requires coordinate-sorted input and refuses paired reads, whose mates would count their overlap twice. Merge mates before alignment, which `cmuts align` does for paired-end input.
- Several alignment files given to one run are read as one merged alignment, where v1 wrote one group per input file. Replicates merge the same way.
- `cmuts sub`, `cmuts div`, and `cmuts norm` read and write whole HDF5 files, where v1's `--experiment` named datasets inside one counts file.
- Each output file holds one flat layout, with `program` and `version` attributes recording what wrote it, in place of the per-file and per-experiment groups and the `meta` group.
- `--min-depth` masks positions below a coverage threshold, absorbing v1's `--blank-cutoff`.
- `cmuts hmm` streams its input and writes no `.cmix` or `.cmfa` index files beside it.

### Removed

- The `sm-dms` and `sm-shape` normalization schemes, and per-reference normalization. `cmuts norm` offers `ubr` and `outlier`, and takes one norm over every input given to a run; to normalize experiments separately, run it once per experiment.
- Termination (RT stop) counting. It may return in a later release.
- The deletion-spreading modes (`--uniform-spread`, `--no-spread`, `--disable-ambiguous`) and `--collapse`, which the pair-HMM subsumes.
- `--max-indel-length`, which was found to hurt performance.
- The filters without a v2 counterpart: `--max-hamming`, `--secondary`, `--downsample`, `--ignore-bases`, and the `--blank-5p` and `--blank-3p` masking.
- The `modification-spectra` dataset, the tokenized sequences under `meta`, and the raw `probability` and `pairwise-coverage` counts.
- Demultiplexing and trimming in `cmuts align`, and the ultraplex, cutadapt, and bowtie2 dependencies.
