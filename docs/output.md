# Outputs

cmuts writes HDF5 files with a consistent schema, containing output reactivity values, statistical errors, and coverage statistics.

## Output Datasets

<!-- BEGIN GENERATED cmuts-hmm FIELDS -->
<div class="field" markdown>

### `coverage`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `zero`

The number of reads in which this base was present, weighted by PHRED scores.

</div>

<div class="field" markdown>

### `reactivity`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

The mutation rate at this base, weighted by PHRED scores and in accordance with the HMM parameters.

</div>

<div class="field" markdown>

### `error`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

Standard error of the reactivity values. Purely the statistical error introduced by finite read depths; does not account for experimental or systemic errors.

</div>

<div class="field" markdown>

### `reads/lengths`

**Shape** `(n, 2l)` · **Type** `uint64` · **Fill** `zero`

The number of reads passing all filters, binned by length.

</div>

<div class="field" markdown>

### `reads/counted`

**Shape** `(n,)` · **Type** `uint64` · **Fill** `zero`

The number of reads passing all filters.

</div>

<div class="field" markdown>

### `reads/rejected`

**Shape** `(n,)` · **Type** `uint64` · **Fill** `zero`

The number of reads rejected by at least one filter.

</div>

<div class="field" markdown>

### `reads/unmapped`

**Shape** `()` · **Type** `uint64` · **Fill** `zero`

The number of reads not aligned to any reference.

</div>
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
