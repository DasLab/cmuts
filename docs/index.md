# Introduction

`cmuts` is a software suite for processing MaP-seq datasets, as well as its cousins RING-MaP and MOHCA-seq, or more generally any experiment where the readouts are mutations to a known reference sequence.

## Installing

Each [release](https://github.com/DasLab/cmuts/releases) carries static `cmuts` binaries for Linux (x86_64, aarch64) and macOS (arm64). These are self-contained and have no dependencies. Alternatively, you can [build from source](from-source.md), which may provide marginal speedups.

The bundled `cmuts-align` helper requires

- [`minimap2`](https://github.com/lh3/minimap2)
- [`samtools`](https://github.com/samtools/samtools)
- [`fastp`](https://github.com/OpenGene/fastp) for paired-end input.

The bundled `cmuts-plot` helper requires a Python 3 installation with the `h5py`, `numpy`, and `plotly` packages.

::::{tab} macOS

```sh
brew install minimap2 samtools fastp
python3 -m pip install h5py numpy plotly
```

Install [Homebrew](https://brew.sh) first if you don't have it.
::::

::::{tab} Debian & Ubuntu

```sh
apt install minimap2 samtools fastp
python3 -m pip install h5py numpy plotly
```

Or load the appropriate modules for your cluster (examples [here](clusters.md)).
::::

## Next Steps

- Learn about the [basics of `cmuts`](basics.md)

```{toctree}
:hidden:

self
from-source
basics
format
```

```{toctree}
:hidden:
:caption: Programs

cmuts-align
cmuts-hmm
cmuts-sub
cmuts-div
cmuts-norm
cmuts-csv
cmuts-plot
cmuts-gen
```

```{toctree}
:hidden:
:caption: Reference

clusters
```
