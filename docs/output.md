# Outputs

cmuts writes HDF5 files with a consistent schema, containing output reactivity values, statistical errors, and coverage statistics.

## Output Datasets

<!-- BEGIN GENERATED cmuts-hmm FIELDS -->
`coverage` — `(n, l)`, float32, fill zero.
The number of reads in which this base was present, weighted by PHRED scores.

`reactivity` — `(n, l)`, float32, fill NaN.
The mutation rate at this base, weighted by PHRED scores and in accordance with the HMM parameters.

`error` — `(n, l)`, float32, fill NaN.
Standard error of the reactivity values. Purely the statistical error introduced by finite read depths; does not account for experimental or systemic errors.

`reads/lengths` — `(n, 2l)`, uint64, fill zero.
The number of reads passing all filters, binned by length.

`reads/counted` — `(n,)`, uint64, fill zero.
The number of reads passing all filters.

`reads/rejected` — `(n,)`, uint64, fill zero.
The number of reads rejected by at least one filter.

`reads/unmapped` — `()`, uint64, fill zero.
The number of reads not aligned to any reference.
<!-- END GENERATED cmuts-hmm FIELDS -->

## Mixed-Length Libraries

For the datasets whose size depends on the reference lengths, the size is always chosen to be the length of the longest reference. For libraries with mixed-length sequences, the rows in the `reactivity`, `error`, and `coverage` datasets are padded with NaN beyond the length of that row's reference.

!!! warning
    NaNs also appear at low-depth bases in the reactivity and error datasets.

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
