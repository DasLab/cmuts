# cmuts

Fast, multithreaded pair-HMM counting of MaP-seq mutations.

## Dependencies

Building and running the C-based programs requires

- A C11 compiler (GCC 6 or newer, or Clang 9 or newer) and `make`,
- [htslib](https://github.com/samtools/htslib) 1.12 or newer,
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) 1.10 or newer.

The bundled alignment helper also requires

- [minimap2](https://github.com/lh3/minimap2),
- [samtools](https://github.com/samtools/samtools).

On macOS with [Homebrew](https://brew.sh):

```sh
brew install htslib hdf5 minimap2 samtools
```

On Debian and Ubuntu:

```sh
apt install libhts-dev libhdf5-dev minimap2 samtools
```

## Installing

```sh
git clone https://github.com/DasLab/cmuts
cd cmuts
make install
```

This installs `cmuts-align`, `cmuts-hmm`, `cmuts-sub`, `cmuts-div`, `cmuts-norm` and `cmuts-gen` to `~/.local/bin`. Optionally pass `BINDIR` to `make` to change the install location.

## Usage

Compute reactivity rates via the pair-HMM on a specific experiment:

```sh
cmuts-hmm -f references.fasta -o treated.h5 treated.bam
```

Subtract rates computed from multiple experiments:

```sh
cmuts-sub -o reactivity.h5 treated.h5 untreated.h5
```

Normalize reactivity rates across experiments:

```sh
cmuts-norm -o apo-normalized.h5 -o holo-normalized.h5 apo.h5 holo.h5
```

## Documentation

See the [docs](https://daslab.stanford.edu/cmuts) for more details.

## Development

See `CONTRIBUTING.md` for information on development builds and running tests.
