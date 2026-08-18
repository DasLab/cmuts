# Building cmuts

## Dependencies

In order to build and run cmuts from source, you need

- A C11 compiler (GCC 6 or newer, or Clang 9 or newer) and `make`,
- [htslib](https://github.com/samtools/htslib) 1.12 or newer,
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) 1.10 or newer.

::::{tab} macOS

```sh
brew install htslib hdf5
```

Install [Homebrew](https://brew.sh) first if you don't have it.
::::

::::{tab} Debian & Ubuntu

```sh
apt install libhts-dev libhdf5-dev
```

Or load the appropriate modules for your cluster (examples [here](clusters.md)).
::::

## Installation

```sh
git clone https://github.com/DasLab/cmuts
cd cmuts
make install
```

This installs the `cmuts` binary and the `cmuts-align` helper script to `~/.local/bin`. Pass `BINDIR` to `make` to install them elsewhere.
