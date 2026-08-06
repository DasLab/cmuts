"""Reduce an output file to the four numbers the tests compare.

Kept and rejected reads, references written, and unmapped reads, in that
order. Nothing here knows what the processing step computes, so these tests
stay valid however it changes.
"""

import sys

import h5py
import numpy as np

with h5py.File(sys.argv[1], "r") as f:
    reads = f["reads"][:]
    rejected = f["reads_filtered"][:]

    print(
        int(np.nansum(reads, dtype=np.float64)),
        int(np.nansum(rejected, dtype=np.float64)),
        int((~np.isnan(reads)).sum()),
        int(f.attrs["reads_unmapped"]),
    )
