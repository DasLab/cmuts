# cmuts

cmuts counts the mutations a chemical probe leaves in MaP-seq reads, and writes a reactivity profile for every reference in an HDF5 file.

It takes coordinate-sorted alignments and the reference sequences they were aligned to:

```sh
cmuts -f references.fasta -o results.h5 alignments.bam
```

Several alignment files are read as one, so a run split across lanes needs no merging first:

```sh
cmuts -f references.fasta -o results.h5 lane1.bam lane2.bam lane3.bam
```

BAM, SAM and CRAM are all read, and no index is needed.

## Installing

cmuts is built against two libraries, and needs a C11 compiler and `make`:

- [htslib](https://github.com/samtools/htslib) 1.12 or newer
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) 1.10 or newer

GCC 6 or newer, or Clang 9 or newer, will compile it.

On macOS with [Homebrew](https://brew.sh), which also installs the Command Line Tools that provide the compiler:

```sh
brew install htslib hdf5
```

On Debian and Ubuntu:

```sh
apt install libhts-dev libhdf5-dev
```

Then:

```sh
git clone https://github.com/hmblair/cmuts
cd cmuts
make install
```

This installs `cmuts` and `cmuts-sub` to `~/.local/bin`. Pass `BINDIR` to `make` to install them elsewhere.

## Where to go next

- [Running cmuts](running.md) — what the options do to a run
- [The output file](output.md) — what a result holds and how to read it
- [Subtracting a background](subtracting.md) — `cmuts-sub`
