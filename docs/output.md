# The output file

A run writes one HDF5 file holding every reference in the FASTA.

## Which row is which reference

Rows are in FASTA order: row *i* holds the *i*th reference in the file that was passed to `--fasta`. The names are not written into the output, so keep the FASTA with the results — it is what identifies the rows, and it holds the sequences besides.

```python
import h5py

def names(fasta):
    return [line[1:].split()[0] for line in open(fasta) if line.startswith(">")]

with h5py.File("results.h5") as f:
    reactivity = f["reactivity"][:]

profile = dict(zip(names("references.fasta"), reactivity))
```

## What each dataset holds

`coverage`, `reactivity` and `error` hold one value per reference position. `reactivity` is the modifications counted at a position over the evidence for them, so it lies between zero and one, and `error` is its standard error. The signal-to-noise at a position is `reactivity / error`.

`reads/lengths` is a histogram of the stored length of every counted read, one row per reference. It is indexed by length rather than by position: column *i* holds the reads of length *i + 1*, and every row is the same width, so a column means the same length in every row. A read longer than the last bin is counted in no bin, which leaves a row summing to fewer reads than `reads/counted`.

`reads/counted` and `reads/rejected` hold one number per reference, and `reads/unmapped` one for the run.

## Missing values

Every row runs to the width the longest reference needs, so most rows are wider than their own reference. **A rate is NaN for two different reasons**: the position is past the end of that row's reference, or it failed `--min-depth`. The coverage tells them apart, being NaN outside a reference and a number within one.

A count is never NaN. A count of zero at a position means the position was reached and nothing was counted there, and a reference no read reached holds a row of zeros rather than a row of NaN.

## Reading a profile

```python
import h5py
import numpy as np

with h5py.File("results.h5") as f:
    reactivity = f["reactivity"][:]
    coverage = f["coverage"][:]

# The bases of one reference, without the padding past its end.
row = reactivity[0][~np.isnan(coverage[0])]
```
