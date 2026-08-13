# Cluster notes

An incomplete list of the clusters cmuts has been built and run on.

## Stanford Sherlock

Running need htslib and HDF5:

```sh
ml load hdf5/1.14.4
ml load biology samtools/1.16.1
```

In addition, building needs a later GCC:

```sh
ml load gcc/12.4.0
```

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
