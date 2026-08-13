# Introduction

cmuts is a software suite for processing MaP-seq datasets, as well as its cousins RING-MaP and MOHCA-seq, or more generally any experiment where the readouts are mutations to a known reference sequence.

## Dependencies

In order to build and run cmuts, you need

- A C11 compiler (GCC 6 or newer, or Clang 9 or newer) and `make`
- [htslib](https://github.com/samtools/htslib) 1.12 or newer
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) 1.10 or newer

=== "macOS"

    ```sh
    brew install htslib hdf5
    ```

    Install [Homebrew](https://brew.sh) first if you don't have it.

=== "Debian & Ubuntu"

    ```sh
    apt install libhts-dev libhdf5-dev
    ```

    Or load the appropriate modules for your cluster (examples [here](clusters.md)).

## Installation

```sh
git clone https://github.com/hmblair/cmuts
cd cmuts
make install
```

This installs all cmuts binaries to `~/.local/bin`. Pass `BINDIR` to `make` to install them elsewhere.

## Next Steps

- Learn about the [basics of cmuts](basics.md)
