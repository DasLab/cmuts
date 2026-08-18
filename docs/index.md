# Introduction

cmuts is a software suite for processing MaP-seq datasets, as well as its cousins RING-MaP and MOHCA-seq, or more generally any experiment where the readouts are mutations to a known reference sequence.

## Installing

Each [release](https://github.com/DasLab/cmuts/releases) carries static `cmuts` binaries for Linux (x86_64, aarch64) and macOS (arm64). These are self-contained and have no dependencies. The bundled `cmuts-align` helper requires

- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)
- [fastp](https://github.com/OpenGene/fastp) for paired-end input.

::::{tab} macOS

```sh
brew install minimap2 samtools fastp
```

Install [Homebrew](https://brew.sh) first if you don't have it.
::::

::::{tab} Debian & Ubuntu

```sh
apt install minimap2 samtools fastp
```

Or load the appropriate modules for your cluster (examples [here](clusters.md)).
::::

Alternatively, you can [build from source](from-source.md), which may provide marginal speedups.

## Next Steps

- Learn about the [basics of cmuts](basics.md)

```{toctree}
:hidden:

self
from-source
basics
output
```

```{toctree}
:hidden:
:caption: Programs

cmuts-align
cmuts-hmm
cmuts-sub
cmuts-div
cmuts-norm
cmuts-gen
```

```{toctree}
:hidden:
:caption: Reference

counting
clusters
```
