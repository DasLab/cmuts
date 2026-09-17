# cmuts

Fast, multithreaded pair-HMM counting of MaP-seq mutations.

## Web Server

The [web server](https://huggingface.co/spaces/DasLab/cmuts) runs `cmuts` on Hugging Face, bypassing the need for a local install. It is suitable for jobs of a few references and a few thousand reads.

## Installing

Each [release](https://github.com/DasLab/cmuts/releases) carries static `cmuts` binaries for Linux (x86_64, aarch64) and macOS (arm64). These bundle htslib and HDF5, so no libraries need to be installed to run them. Alternatively, you can [build from source](https://daslab.stanford.edu/cmuts/from-source), which may provide marginal speedups.

`cmuts align` calls separate programs, which must be on the `PATH`.

- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)
- [vsearch](https://github.com/torognes/vsearch) for paired-end input.

On macOS with [Homebrew](https://brew.sh):

```sh
brew install minimap2 samtools vsearch
```

On Debian and Ubuntu:

```sh
apt install minimap2 samtools vsearch
```

`cmuts plot` requires a Python 3 installation with the `h5py`, `numpy`, and `plotly` packages.

```sh
python3 -m pip install h5py numpy plotly
```

## Usage

Compute reactivity rates via the pair HMM on a specific experiment:

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

Generate an interactive report of the results:

```sh
cmuts plot \
    apo-normalized.h5 --label "Apo" \
    holo-normalized.h5 --label "Holo"
```

## Documentation

See the [docs](https://daslab.stanford.edu/cmuts) for more details on each program, the HDF5 outputs, and special use cases.

## Development

See [`CONTRIBUTING.md`](https://github.com/DasLab/cmuts?tab=contributing-ov-file) for information on development builds and running tests.
