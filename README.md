# cmuts

Fast, multithreaded pair-HMM counting of MaP-seq mutations.

## Dependencies

- A C11 compiler and `make`
- [htslib](https://github.com/samtools/htslib)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)

On macOS with Homebrew:

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

This installs `cmuts` to `~/.local/bin`. Optionally set `BINDIR` to change the install location.
