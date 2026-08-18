# cmuts

Fast, multithreaded pair-HMM counting of MaP-seq mutations.

## Installing

Each [release](https://github.com/DasLab/cmuts/releases) carries static `cmuts` binaries for Linux (x86_64, aarch64) and macOS (arm64). These are self-contained and have no dependencies. The bundled `cmuts-align` helper requires

- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)
- [fastp](https://github.com/OpenGene/fastp) for paired-end input.

On macOS with [Homebrew](https://brew.sh):

```sh
brew install minimap2 samtools fastp
```

On Debian and Ubuntu:

```sh
apt install minimap2 samtools fastp
```

Alternatively, you can [build from source](https://daslab.stanford.edu/cmuts/from-source), which may provide marginal speedups.

## Usage

Compute reactivity rates via the pair-HMM on a specific experiment:

```sh
cmuts hmm -f references.fasta -o treated.h5 treated.bam
```

Subtract rates computed from multiple experiments:

```sh
cmuts sub -o reactivity.h5 treated.h5 untreated.h5
```

Normalize reactivity rates across experiments:

```sh
cmuts norm -o apo-normalized.h5 -o holo-normalized.h5 apo.h5 holo.h5
```

## Documentation

See the [docs](https://daslab.stanford.edu/cmuts) for more details.

## Development

See `CONTRIBUTING.md` for information on development builds and running tests.
