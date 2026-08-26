# File Format

Every HDF5 file `cmuts` reads and writes follows a consistent format, containing or expecting a subset of the datasets below. See an individual program's page for the datasets it reads or writes and how it produces them.

The size of each dataset depends on the number of references `n` and the length of the longest reference `l`. Ragged libraries contain left-aligned data, with the fill value where the reference did not have any bases.

The type indicated only applies when the corresponding dataset is written, as numeric types are automatically converted during loading.

## Datasets

<!-- BEGIN GENERATED cmuts FIELDS -->
{.field}
### `mismatches/rate`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

The rate of mismatches at each base, over the coverage.

{.field}
### `mismatches/error`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

Binomial standard error of the mismatch rate.

{.field}
### `insertions/rate`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

The rate of insertions opened after each base, over the reads pairing the base and continuing 5'.

{.field}
### `insertions/error`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

Binomial standard error of the insertion rate.

{.field}
### `deletions/rate`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

The rate of deletion runs ending at each base, over the reads pairing the base 3' of it and continuing 5'.

{.field}
### `deletions/error`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `NaN`

Binomial standard error of the deletion rate.

{.field}
### `norm`

**Shape** `()` · **Type** `float32` · **Fill** `NaN`

The norm every rate in this file was divided by.

{.field}
### `coverage`

**Shape** `(n, l)` · **Type** `float32` · **Fill** `0`

The number of reads in which each base was present.

{.field}
### `sequence`

**Shape** `(n, l)` · **Type** `int8` · **Fill** `-1`

The reference sequence: 0 for A, 1 for C, 2 for G, 3 for T, and -1 for any other base and for every column past the reference's end.

{.field}
### `reads/counted`

**Shape** `(n,)` · **Type** `uint64` · **Fill** `0`

The number of reads contributing to the rates.

{.field}
### `reads/lengths`

**Shape** `(n, 2l)` · **Type** `uint64` · **Fill** `0`

The number of reads contributing to the rates, binned by length.

{.field}
### `reads/rejected`

**Shape** `(n,)` · **Type** `uint64` · **Fill** `0`

The number of reads aligned to the reference but not contributing to the rates.

{.field}
### `reads/unmapped`

**Shape** `()` · **Type** `uint64` · **Fill** `0`

The number of reads not aligned to any reference.

{.field}
### `pairwise/correlation`

**Shape** `(n, l, l)` · **Type** `float32` · **Fill** `NaN`

The Pearson correlation of mutations between this pair of bases.

{.field}
### `pairwise/conditional`

**Shape** `(n, l, l)` · **Type** `float32` · **Fill** `NaN`

The probability that the base on the first axis was mutated in a read, given that the base on the second axis was.

{.field}
### `pairwise/coverage`

**Shape** `(n, l, l)` · **Type** `float32` · **Fill** `0`

The number of reads in which each pair of bases was present.
<!-- END GENERATED cmuts FIELDS -->

## File Attributes

Every file records metadata on the root group. These are attributes and not datasets, so `h5py` reads them from `f.attrs`.

<!-- BEGIN GENERATED cmuts ATTRIBUTES -->
| Attribute | Description |
| --- | --- |
| `program` | The name of the program that produced this file. |
| `version` | The version of `cmuts` that produced this file. |
<!-- END GENERATED cmuts ATTRIBUTES -->

## Reading With Python

```python
import h5py

# Reading attributes

with h5py.File("output.h5") as f:
    print(f.attrs["program"], f.attrs["version"])

# Reading datasets

with h5py.File("output.h5") as f:
    reactivity = f["reactivity"][:]
    reads = f["reads/counted"][:]
```
