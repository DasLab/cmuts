# cmuts

Fast, multithreaded pair-HMM counting of MaP-seq mutations.

## Dependencies

- A C11 compiler (GCC 6 or newer, or Clang 9 or newer) and `make`
- [htslib](https://github.com/samtools/htslib) 1.12 or newer
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) 1.10 or newer

On macOS with [Homebrew](https://brew.sh):

```sh
brew install htslib hdf5
```

On Debian and Ubuntu:

```sh
apt install libhts-dev libhdf5-dev
```

## Installing

```sh
git clone https://github.com/hmblair/cmuts
cd cmuts
make install
```

This installs `cmuts-hmm`, `cmuts-sub`, `cmuts-div`, `cmuts-norm` and `cmuts-gen` to `~/.local/bin`. Optionally pass `BINDIR` to `make` to change the install location.

## Development

See `CONTRIBUTING.md` for information on development builds and running tests.
