# Outputs

cmuts writes HDF5 files with a consistent schema, containing output reactivity values, statistical errors, and coverage statistics.

## File Attributes

Every file records what produced it, on the root group. These are attributes and not datasets, so `h5py` reads them from `f.attrs`.

<!-- BEGIN GENERATED cmuts-hmm ATTRIBUTES -->
| Attribute | Description |
| --- | --- |
| `program` | The name of the program that produced this file. |
| `version` | The version of cmuts that produced this file. |
<!-- END GENERATED cmuts-hmm ATTRIBUTES -->

```python
with h5py.File("output.h5") as f:
    print(f.attrs["program"], f.attrs["version"])
```

## Output Datasets

<!-- BEGIN GENERATED cmuts-hmm FIELDS -->
{.field}
### `coverage`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `zero`

The number of reads in which this base was present, weighted by PHRED scores.

{.field}
### `reactivity`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

The mutation rate at this base, weighted by PHRED scores and in accordance with the HMM parameters.

{.field}
### `error`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

Standard error of the reactivity values. Purely the statistical error introduced by finite read depths; does not account for experimental or systemic errors.

{.field}
### `reads/lengths`

**Shape** `(n, 2l)` · **Type** `uint64` · **Fill** `zero`

The number of reads passing all filters, binned by length.

{.field}
### `reads/counted`

**Shape** `(n,)` · **Type** `uint64` · **Fill** `zero`

The number of reads passing all filters.

{.field}
### `reads/rejected`

**Shape** `(n,)` · **Type** `uint64` · **Fill** `zero`

The number of reads that contributed nothing: rejected by at least one filter, or carrying something the pair HMM's rates give no alignment of.

{.field}
### `reads/unmapped`

**Shape** `()` · **Type** `uint64` · **Fill** `zero`

The number of reads not aligned to any reference.

{.field}
### `pairwise/correlation`

**Shape** `(n, l, l)` · **Type** `float32` · **Fill** `NaN` · _written only when asked for_

The correlation between two positions being modified in the same read, as the Pearson coefficient of the two binary variables. NaN where the reads are too few, and where either position is modified in all of them or in none. The diagonal is a position against itself, which falls short of one by however much of its variance the base calls leave unsettled; divide a correlation by the square root of the two diagonals to take that out.

{.field}
### `pairwise/coverage`

**Shape** `(n, l, l)` · **Type** `float32` · **Fill** `zero` · _written only when asked for_

The evidence behind each correlation: the reads reaching both positions.

{.field}
### `norm`

**Shape** `()` · **Type** `float32` · **Fill** `NaN` · _written only when asked for_

The scale every rate in this file was divided by.
<!-- END GENERATED cmuts-hmm FIELDS -->

## Mixed-Length Libraries

For the datasets whose size depends on the reference lengths, the size is always chosen to be the length of the longest reference. For libraries with mixed-length sequences, the rows in the `reactivity`, `error`, and `coverage` datasets are padded with NaN beyond the length of that row's reference.

```{warning}
NaNs also appear at low-depth bases in the reactivity and error datasets.
```

## Reading Outputs

```python
import h5py
import numpy as np

with h5py.File("output.h5") as f:
    reactivity = f["reactivity"][:]
    error = f["error"][:]
```

### Signal-To-Noise

```python
snr = reactivity / error                     # per-base
snr = np.mean(reactivity / error, axis=-1)   # per-reference
```
