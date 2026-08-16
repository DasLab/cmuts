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

Which datasets a file holds depends on what wrote it, and so does what the numbers in
them mean. Each program's page lists the datasets it writes and describes them as they
stand in its own output.

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
