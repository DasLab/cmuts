# Changelog

## [1.5.0] - 2026-07-29

### Added

- `cmuts normalize --norm sm-dms`: a ShapeMapper2-style per-nucleotide DMS normalization scheme. Each base (A, C, G, U) is scaled independently by the 75th percentile of that base's reactivities pooled across the whole library (matching ShapeMapper2's single `normalize_profiles` pass), with each position divided by its own base's factor. Sequence-aware schemes now receive the per-position base identity from the required `--fasta` when the counts file carries no tokenized sequence
- `cmuts normalize --norm sm-shape`: a ShapeMapper2-style SHAPE boxplot normalization scheme, reproducing ShapeMapper2's `find_boxplot_factor` exactly (drop outliers above `max(1.5*IQR, 90th-percentile value)`, then divide by the mean of the top 10% of the remainder; the 95th percentile is used below 100 values). With `sm-dms` (DMS) it lets cmuts reproduce ShapeMapper2's normalization for both probes
- `cmuts core` now reads SAM alignment files (plain text or bgzipped), alongside BAM and CRAM. SAM records decode to the same internal alignment as BAM via htslib's `sam_parse1`, reusing the existing per-reference index and seek machinery (BGZF reads uncompressed files transparently). Format is detected from the file's leading bytes rather than the extension. SAM input must be coordinate-sorted, as with BAM and CRAM. The test matrix is now parametrized over all three formats
- A `benchmarks/` suite comparing cmuts against rf-count (RNAFramework) and shapemapper2, built on a shared package (`benchmarks/external/`, one module per tool) that is the single place that runs each tool's full reactivity pipeline. Its entry point, `reactivity(tool, count, norm, inputs, condition, out_h5)`, takes an `Inputs` struct plus the count- and normalize-stage parameter structs (`CountParams`/`NormParams`), adapts each input alignment to the format the tool needs (cmuts reads bam/cram/sam natively; rf-count converts to BAM; shapemapper2 sorts/indexes to an indexed BAM), counts, normalizes, and writes a unified HDF5 with per-reference reactivity and read depth (the max post-filter per-position coverage, matching cmuts' own measure) — applying no coverage floor of its own, so consumers filter as they choose. Three scripts: `profile` (wall-clock time and peak memory of each tool's full pipeline across query/reference/length sweeps, on its own synthetic data, for BAM/CRAM/SAM), `tool-vs-tool` (per-nucleotide agreement of cmuts-tuned-to-each-tool against that tool — the control that isolates deletion spreading), and `tool-vs-structure` (AUC/Pearson/Spearman of reactivity against a base-pairing ground truth read from the deposited mmCIF hydrogen-bond annotations). The two scoring scripts build the reactivity they need directly from the modified/untreated alignments via `external.reactivity`. Documented on the new Benchmarks docs page
- A per-base sequence view in the `cmuts plot` report: each base of the reference is drawn in a monospace strip colored from white (low) to purple (high) by its reactivity, with masked positions muted. It appears under the summary statistics for single-reference experiments and under the per-reference statistics — tracking the selected reference — for multi-reference experiments (`--all`). Driven by the new `meta/sequence` dataset (decoded from base tokens to A/C/G/U client-side)

### Changed

- `cmuts generate` now generates synthetic test data in parallel across references with OpenMP, giving a ~6.5x speedup on a 15-core machine. Each reference is independent — it fills its own slice of the in-memory HDF5 buffer and draws from its own generator seeded `base_seed + ref_ix`, so output is reproducible and bit-identical regardless of thread count. The per-chunk HDF5 write stays serial (no HDF5 thread-safety or MPI-collective I/O required), and the FASTA is written in reference-index order so it stays aligned with the SAM `@SQ` header and the expected-value layout that `cmuts core` keys on. The executable build targets now define `HAVE_OPENMP` (previously they linked OpenMP but left the macro undefined, so the generator's parallel regions compiled out)
- `cmuts normalize` now specifies inputs with a single repeatable `--experiment NAME mod=... [nomod=...]` flag (breaking). This replaces the separate single-group (`--mod`/`--nomod`/`--group`) and multi-group (`--groups` TOML) interfaces — one experiment or many use the same syntax, and the TOML format is removed. The experiment name is the output HDF5 group. In the Python API, `Group`/`GroupResult` are renamed to `Experiment`/`ExperimentResult`
- Normalization granularity is now controlled by two composable flags, `--per-experiment-norm` and `--per-reference-norm`, replacing `--independent-norm` (breaking). By default a single factor is shared across all experiments and references (maximally comparable); `--per-experiment-norm` normalizes each experiment on its own, and `--per-reference-norm` normalizes each reference on its own (per-transcript, like rf-norm). Both axes compose with every `--norm` scheme — e.g. `outlier` with `--per-reference-norm` reproduces rf-norm's per-transcript 2-8%, without it the library-wide variant
- Plotting is decoupled from normalization (breaking). `cmuts normalize` now only computes: it writes the reactivity HDF5 (a complete record of the plottable data) and prints statistics, but no longer generates figures. Run `cmuts plot <reactivity.h5>` to render them — `cmuts normalize` followed by `cmuts plot` reproduces the previous behavior
- `cmuts plot` now generates a single self-contained, interactive HTML report from a reactivity HDF5, replacing the directory of static image files (breaking). The report is driven entirely by the HDF5 dataset names (the contract between `cmuts normalize` and the report — no code is shared), with collapsible per-experiment sections: summary statistics, per-base mutation heatmap, coverage, terminations, SNR-vs-depth, and (for multi-reference experiments) a reference-clipped reactivity heatmap and read-distribution plots. `cmuts plot --all` additionally embeds a per-sequence dropdown (that reference's profile, coverage, terminations, mutual information, and correlation) plus two cross-experiment comparison sections — an arbitrary-overlay multi-plot and a two-reference difference plot. The interface is now `cmuts plot <file.h5> [--all] [--group NAME] [-o report.html]`, where `-o` names an HTML file rather than a directory
- Plotly is now the sole plotting backend and a core dependency; matplotlib is removed as a dependency and the `hf` optional-dependency extra (which previously carried Plotly) is dropped (breaking). Every figure in the report is a Plotly plot, and static images come from Plotly's modebar "Download plot as PNG" button at a fixed, window-independent size — so there is no headless-browser (Chrome/kaleido) dependency
- Shared, experiment-independent datasets now live under a top-level `meta` group: `meta/sequence` (per-reference base tokens) and `meta/name` (per-reference names in FASTA order), replacing the previous top-level `sequence`/`names` datasets (breaking). `cmuts normalize` now always writes the reference sequence (previously only when the active `--norm` scheme consumed it), so the report can color bases by reactivity. `meta` is reserved and rejected as an experiment name; readers fall back to the old top-level datasets for files written before the move
- Plotly builder API (`cmuts.visualize`) reworked to take raw arrays and match the report's needs (breaking): `plot_snr_scaling` takes the raw SNR curve arrays instead of an `SNRCurves` object; `plot_termination` takes a single reference's raw 1-D termination counts and normalizes internally; and `plot_multiple_examples` is renamed `plot_reactivity_heatmap`
- Removed `cmuts normalize --norm-cutoff` and `--norm-percentile` (breaking). They were parsed but never applied; the UBR scheme uses a 500-read coverage floor and the 90th percentile (the documented defaults)
- Removed the public `NormScheme` enum (breaking). Normalization schemes are now a registry from which the `--norm` choices are derived, so adding a scheme is a single registered class

### Fixed

- `configure --mpi` failed to compile: including `<mpi.h>` pulled in OpenMPI's deprecated C++ bindings, whose `namespace MPI` collided with cmuts' own `MPI` namespace (redeclaration of `MPI::Init_thread`). The MPI header now defines `OMPI_SKIP_MPICXX`/`MPICH_SKIP_MPICXX` before the include to suppress those bindings, and drops the unneeded `extern "C"` wrapper (`mpi.h` supplies its own linkage guards)
- `cmuts normalize` could be OOM-killed on reference-heavy datasets: the SNR-vs-read-depth plot materialized a `(grid × references × positions)` array (tens of GB at 10,000 grid points across thousands of references). The SNR curves are now computed once in the data layer with a memory-bounded reduction and a 1,000-point grid, and saved into the HDF5 for the plotter to render
- Shell completions for `normalize`, `plot`, and `visualize` had drifted from the actual CLIs (they offered a nonexistent `--norm-independent`, and omitted `--groups` and several `cmuts visualize` flags). Completions are now generated from the argparse parsers via shtab, so they stay in sync
- The reactivity HDF5 now stores per-reference, per-position coverage and per-position termination densities (previously reconstructed or zeroed on load), so `cmuts plot` reproduces the coverage and termination figures faithfully — including a single reference's coverage in the per-sequence view
- `cmuts core` now handles alignment or FASTA files in read-only directories. The per-file `.cmix` and `.cmfa` index files are normally written beside their source; when that directory is not writable, the index is redirected into the system temp directory under a name derived from the absolute source path (unique per source and reused across runs). Previously a read-only alignment directory caused the file to be silently dropped (cmuts exited 0 with empty output) and a read-only FASTA directory caused an infinite hang while parsing the missing index. Index writes that genuinely fail now report a clear error, and a missing/truncated FASTA index is detected instead of looping
- `configure` no longer leaves `bin/` missing the binaries or dispatcher symlinks. It previously wiped `bin/` unconditionally at the start (so a failed build left no working install) and created the `cmuts`/`cmuts-generate`/... symlinks only after the `--build-only` early-exit (so `--build-only` produced binaries but no dispatcher commands). The `bin/` wipe is now gated on `--clean`, and the symlinks are created before the `--build-only` gate with idempotent `ln -sf`, so `bin/` is always fully assembled
- The `cmuts plot` reactivity-difference plot no longer draws spurious values across masked positions. The difference is computed in the browser as `a - b`, and masked positions (embedded as JSON `null`) were coercing to `0` in the subtraction (`null - x === -x`); either operand being masked now yields a gap, so the line and its error band break there like every other profile

### Removed

- The matplotlib plotting backend (`cmuts.visualize.plotting`) and its figure functions. Plotting is Plotly-only; `cmuts-plot` now delegates to `cmuts.report:main`
- The Plotly builders `plot_examples`, `plot_profile`, `plot_profiles`, `plot_error`, and `plot_variance`, and the `sequence` parameter added to `plot_profile`/`plot_profiles`/`plot_examples` in 1.4.8. The report renders profiles client-side (a single shared line-plus-error-band function), and the per-position nucleotide axis is superseded by the colored sequence view
- The unused `cmuts.requires_sequence` helper (its only caller — the gate that attached the sequence only for sequence-aware norm schemes — was removed now that the sequence is always stored). `Scheme.needs_sequence`, which decides whether a norm scheme consumes base identity, is unchanged

## [1.4.10] - 2026-06-22

### Added

- `cmuts core --token-map` flag for assigning custom non-negative integer tokens to A, C, G, U in the `--tokenize` output (e.g. `--token-map 4,5,6,7`). Defaults to `0,1,2,3`; T and U share the fourth value

### Changed

- Renamed `cmuts normalize` flags `--clip-low`/`--clip-high` to `--clip-below`/`--clip-above` (breaking)
- Renamed `plot_cumulative_reads` to `plot_reads_per_block` (breaking)
- `cmuts core --downsample` now caps the number of reads that *pass* the quality filters per reference, rather than the number of reads read from the file. This bounds usable depth and avoids over-processing references
- `cmuts -h` now lists the subcommands with one-line descriptions, and `cmuts test -h` shows the test-runner options (previously the test options were shown under the top-level help)

### Fixed

- `cmuts core` now decodes CRAM files written by real aligners (e.g. Ultima `sorter`). The hand-rolled CRAM decoder previously aborted such files with "Invalid read length" or "reference position exceeds reference length": it never decoded multi-symbol Huffman codes (it returned a constant and consumed no bits) and read the per-record fields out of CRAM-spec order. Canonical Huffman decoding is now implemented and the fields are read in spec order. Records that use CRAM features cmuts does not decode (core-coded mate or auxiliary-tag data) now fail cleanly instead of silently mis-decoding
- `cmuts core --tokenize` no longer writes all-zero tokens for reference sequences shorter than four bases (a leftover stream-error flag from header parsing turned the binary FASTA read into a silent no-op)
- `cmuts core --tokenize` now exits non-zero when tokenization fails, such as on an invalid `--token-map`

## [1.4.9] - 2026-05-08

### Added

- `cmuts.visualize.make_defattr` and `cmuts.visualize.chimerax_command` public functions for splitting structure visualization into "write defattr" and "build ChimeraX command" steps — enables server-side workflows where the user runs ChimeraX locally
- `visualize_structure` now accepts `defattr_out` and `run_chimerax` parameters and returns the ChimeraX command string
- `cmuts-visualize` CLI gains `--defattr-out` and `--no-chimerax` flags. With `--no-chimerax`, the script writes the defattr and prints the ChimeraX command instead of launching ChimeraX

## [1.4.8] - 2026-05-08

### Added

- Optional `sequence` parameter on `plot_profile`, `plot_profiles`, and `plot_examples` (Plotly). When provided, the x-axis shows the nucleotide letter at every position with the residue number underneath at every 10th position (and the first/last)

## [1.4.7] - 2026-05-06

### Added

- `cmuts.Group` dataclass and `cmuts.compute_reactivities()` for processing multiple experimental groups in one pass, with normalization pooled across groups so reactivities are directly comparable
- `cmuts.save_groups()` for writing multiple groups to a single HDF5 file
- `cmuts.pooled_norm()` for computing a shared normalization factor across multiple `ProbingData` objects (UBR pools high-coverage reactivities; OUTLIER averages per-reference factors)
- `cmuts normalize --groups <toml>` flag accepting a TOML file with `[[group]]` entries (`name`, `mod`, optional `nomod`)
- `cmuts normalize --independent-norm` flag to disable cross-group pooling
- Multi-group integration tests for the `cmuts normalize` CLI
- Plotly plotting backend (`cmuts.visualize.plotly`) for interactive web displays, including profile, pairwise heatmap, and examples plots
- `--vmin`/`--vmax` clipping flags on `cmuts-visualize`

### Changed

- `cmuts-normalize` script is now a thin shim over `python -m cmuts.normalize`; the implementation lives in `cmuts/normalize/__main__.py`
- Pairwise Plotly heatmaps now enforce a 1:1 aspect ratio with cell outlines, hover text, and constrained layout

### Fixed

- `cmuts-visualize` CLI no longer fails on missing arguments
- HDF5 dataset-not-found errors now report the dataset name and available keys
- Build no longer fails when `setuptools-scm` cannot determine a version (fallback version added)
- GitHub Actions workflows updated to Node 24-ready action versions

## [1.4.6] - 2026-03-17

### Added

- SNR scaling plot for `cmuts normalize`, showing how mean SNR changes with relative sequencing depth
- Optimal read-allocation ("Pareto") curve and infeasible-region hatching for two-condition SNR scaling plots
- Error bands on SNR scaling curves to visualize uncertainty in projected mean SNR

### Fixed

- Plot outputs now save consistently to the requested output directory, including single-profile plots
- Multi-sequence pairwise plots now use distinct filenames instead of overwriting one another
- Plotting now uses a headless matplotlib backend by default, avoiding aborts in non-interactive environments
- `cmuts normalize` test harness now uses the active Python interpreter, improving reliability in virtual environments

### Changed

- `cmuts normalize` now uses signed SNR in reports and visualizations
- Normalize documentation updated for the signed-SNR workflow and new SNR scaling plot
- Figures generated by `cmuts normalize` are written next to the output HDF5 file instead of relative to the current working directory
- Build and packaging portability improved, including Bioconda recipe updates and htslib/CMake discovery cleanup

## [1.4.3] - 2026-03-12

### Added

- Bash and zsh shell completions for the `cmuts` CLI, with context-aware file type filtering for all subcommands
- `--install-completions` flag for `./configure` to install completions to user directories
- `scripts/generate_completions.py` to generate completions from JSON arg configs
- GitHub workflow to sync cmuts source to HF Space repo on push
- Hugging Face Space link in README and docs

### Fixed

- Memory units label corrected from KB to MB in profile script

### Changed

- Progress updates suppressed in non-TTY mode; only final stats are printed
- Terminal escape sequences suppressed when stdout is not a TTY

## [1.4.2] - 2026-02-24

### Added

- Integration tests for `cmuts normalize` CLI script
- `--sam` flag to `cmuts generate` for explicit SAM output
- ShapeMapper2 benchmarking to profile script

### Fixed

- `cmuts-normalize` script updated to match v1.4.0 Python API:
  - `cmuts.normalize()` renamed to `cmuts.compute_reactivity()`
  - `cmuts.visualize.all()` renamed to `cmuts.visualize.plot_all()`
  - `--raw` and `--norm-outlier` flags replaced with `--norm {ubr,raw,outlier}`
  - `Opts()` call updated to pass `norm` string instead of `raw` bool
- `ProbingData.normalize()` now uses numpy instead of dask operations, fixing compatibility with older dask versions
- Outlier normalization factor now broadcasts correctly against 2D reactivity arrays
- Profile script: removed `--low-mem`, fixed `--cram` flag, use default `max-indel-length`
- Skip `mpirun` for single-thread runs and validate `--threads`
- Filter unaligned reads from SAM before ShapeMapper2 benchmark

### Changed

- `cmuts generate` now requires explicit format selection (`--bam`, `--cram`, or `--sam`)
- Parallelized build in configure script
- CI: removed unused Codecov upload step, fixed PATH expansion and integration test execution

## [1.4.1] - 2026-02-23

### Added

- `--bam` and `--cram` flags to `cmuts generate` for producing sorted BAM/CRAM output directly (requires samtools). Both flags can be used together to produce both formats in one invocation.
- Automated docs deployment to GitHub Pages

### Fixed

- CRAM reader now reads `primary` and `reversed` flags from the BF field instead of hardcoding them
- CLI error messages now show the argument name (e.g., `--length: required.` instead of `: required.`), caused by empty short names being passed to argparse

### Changed

- **Breaking**: `--out-sam` renamed to `-o`/`--out` in `cmuts generate`, since output is no longer necessarily SAM
- Build outputs executables directly to `bin/` instead of requiring a separate install step
- `cmuts-generate-tests` binary renamed to `_cmuts-generate-tests` (internal); users should use `cmuts generate` instead
- Test harness and profile script now use `cmuts generate --bam`/`--cram` instead of manual samtools post-processing

## [1.4.0] - 2025-01-23

### Changed

- **Breaking**: Consolidated normalization CLI flags into single `--norm {ubr,raw,outlier}` argument
  - Replaces `--raw` and `--norm-outlier` flags
  - Default is `ubr` (90th percentile normalization)
- **Breaking**: `Opts` dataclass now uses `norm: str` instead of `raw: bool` and `outlier: bool`

## [1.3.3] - 2025-01-23

### Changed

- **Python normalization extraction**: Moved normalization code from `internal.py` into `normalize/schemes.py`
  - `NormScheme` enum and `get_norm()` function now in `cmuts.normalize.schemes`
  - Reduces `internal.py` from 921 to 831 lines

## [1.3.2] - 2025-01-22

### Changed

- **C++ logging extraction**: Extracted logging code from `utils.hpp/cpp` into dedicated `infra/logging.hpp/cpp` module
  - Logging constants, macros, and functions now in separate header
  - `utils.hpp` includes `logging.hpp` for backward compatibility

## [1.3.1] - 2025-01-22

### Changed

- **Codebase reorganization** for improved maintainability:
  - Python: Extracted formatting/display code into new `cmuts.output` module
  - Python: Renamed `visualize/vis.py` to `visualize/structure.py` for clarity
  - C++: Reorganized into `io/`, `core/`, `infra/`, `app/` subdirectories
  - C++: Split large `common.hpp` into modular headers (`io/common/types.hpp`, `stream.hpp`, `alignment.hpp`, `file.hpp`)
  - Added backward compatibility headers in `compat/`

### Added

- `cmuts.output` module with `stats()`, `title()`, `subtitle()` functions
- `cmuts.subtitle()` function for statistics headers

## [1.3.0] - 2025-01-22

### Added

- CI/CD infrastructure with GitHub Actions workflow
- Automated testing (pytest) and coverage reporting
- Documentation generation in CI pipeline

### Fixed

- CI: Install C++ binaries and fix coverage path
- CI: Add libbz2-dev and liblzma-dev for htscodecs build
- CI: Build htscodecs submodule before cmake
- CI: Update xtensor to 0.27.0 and xtl to 0.8.0

## [1.2.4] - 2025-01-22

### Added

- Configurable logging system via `CMUTS_LOG_LEVEL` and `CMUTS_LOG_STDERR` environment variables
- `TRACE` log level for detailed MPI debugging
- Pairwise SNR calculation for joint probability P(i=1, j=1) in `cmuts normalize --pairwise`

### Fixed

- **Critical**: MPI deadlock in HDF5 collective I/O when processes had asymmetric workloads
- Root cause: `skip` flag was reset after each read, causing processes to participate in different numbers of collective writes

### Changed

- Version now derived from git tags (configure generates VERSION file)

## [1.2.3] - 2025-12-01

### Added

- Pytest integration test suite (`cmuts test`) with edge case, CRAM, fuzzing, and stress tests
- `--verbose` flag to `cmuts core` for debug logging
- Structured logging with levels (DEBUG, INFO, WARN, ERROR) to `.cmuts.log`
- `Alignment::create()` factory method to prevent field order initialization bugs
- Validation macros (`CMUTS_CHECK_BOUNDS`, `CMUTS_CHECK_NOT_NULL`, etc.)
- Linux installation instructions (Ubuntu/Debian, Fedora/RHEL) to documentation
- `cmuts test` documentation page

### Fixed

- **Critical**: CRAM field order bug where `primary` and `reversed` fields were swapped, causing all CRAM reads to be filtered as non-primary
- Linux compatibility in profile script (portable GNU time detection)

### Changed

- `cmuts test` now runs pytest tests by default (use `--quick` for fast tests, `--cram` for CRAM-only)
- Documentation: fixed outdated `cmuts pairwise` reference in pipeline.md

### Removed

- Legacy bash-based test script (`cmuts-test`, `cmuts-test-compare`)

## [1.2.2] - 2025-11-19

### Added

- Fast-path for processing sparse alignments (reads << references)
- `--secondary` flag for processing secondary alignments
- JSON-based argument file for `cmuts core` and `cmuts generate`

### Changed

#### Non-Breaking

- Minor bugfixes
- Reads with missing mapping qualities are skipped, unless `--min-mapq` is set to `0`
- Seek BAM/CRAM only on MPI builds for slightly faster processing
- Replaced file-based mutex with semaphores

## [1.2.1] - 2025-11-14

### Added

- New subprogram `cmuts plot` to plot multiple reactivity profiles, or profiles generated outside of `cmuts normalize`
- Re-instated `--print-every` argument of `cmuts core`, now measured in seconds rather than reads
- Added `--downsample` argument to `cmuts core`, to limit read depths per reference
- Lock for when multiple processes attempt to generate `.cmfa` or `.cmix` files

### Changed

#### Non-Breaking

- Upgraded algorithm for detecting and spreading ambiguous deletions
- `cmuts normalize` now compatible with Python 3.9 and Python 3.16

## [1.2.0] - 2025-11-07

### Added

- Strand filtering via `--no-reverse` and `--only-reverse`
- Hamming distance filtering via `--max-hamming`
- Base filtering via `--ignore-bases`, such as for DMS probing
- Speedup of pairwise profile processing via OMP
- `Low-quality bases` stat added to `cmuts core`
- Length-based blanking for libraries with variable-length references in `cmuts normalize`
- Additional error and Bonferroni-corrected correlation plots in `cmuts normalize`
- Optional `--sig` flag in `cmuts normalize` to specify correlation significance levels
- `--no-insertions` and `--no-deletions` added to `cmuts normalize`
- Python API to `cmuts visualize` (`cmuts.visualize`) and `cmuts normalize` (`cmuts.normalize`)
- `--clean` flag for `configure` script to build from scratch

### Changed

#### Breaking

- `--joint` argument of `cmuts core` renamed `--pairwise`
- Output datasets from `cmuts core` now named `counts-1d` and `counts-2d` for normal and `--pairwise` mode, respectively
- Merged `cmuts pairwise` into `cmuts normalize`, which will auto-process 2D counts if present alongside 1D counts
- `cmuts normalize` now requires a `--fasta` input

#### Non-Breaking

- 1D counts computed even when running `cmuts core` in `--pairwise` mode
- `configure` script will `pip`-install the `cmuts` Python API

### Removed

- Unsused `--subsample` functionality
- `--print-every` argument (now time-based printing)
- Unsused and only partially supported `--low-mem` runtime mode
- Now-merged `cmuts pairwise`
