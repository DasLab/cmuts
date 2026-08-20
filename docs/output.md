# Outputs

cmuts writes HDF5 files with a consistent schema, containing output reactivity values, statistical errors, and coverage statistics.

## File Attributes

Every file records metadata on the root group. These are attributes and not datasets, so `h5py` reads them from `f.attrs`.

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

## Mixed-Length Libraries

The size of length-dependent datasets is always based on the length of the longest reference. For libraries with mixed-length sequences, the positions beyond the length of a row's reference hold the dataset's fill value.

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
