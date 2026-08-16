# Cluster notes

An incomplete list of the clusters cmuts has been built and run on.

## Stanford Sherlock

Running needs htslib and HDF5:

```sh
ml load hdf5/1.12.0
ml load biology samtools/1.16.1
```

In addition, building needs a later GCC:

```sh
ml load gcc/12.4.0
```

The HDF5 modules marked `(m)` in `module avail`, 1.14.4 among them, are built against MPI and bring UCX in with it. UCX installs memory hooks from a constructor, which a sanitized binary crashes on before reaching `main`. 1.12.0 is a serial build and has neither.

## HHMI Janelia

Running needs htslib on the path and its libraries alongside:

```sh
ml load samtools
export LD_LIBRARY_PATH=/misc/local/samtools-1.22.1/lib:$LD_LIBRARY_PATH
```

Building needs htslib where `pkg-config` will find it:

```sh
export PKG_CONFIG_PATH=/misc/local/samtools-1.22.1/lib/pkgconfig:$PKG_CONFIG_PATH
```
