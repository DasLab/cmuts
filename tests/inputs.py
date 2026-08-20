"""Values written into an input file by hand.

What a combination writes depends on the values in its inputs and on nothing
else, so the files it is tested against are written by hand and not counted from
an alignment. outputs.py describes the layout they share.
"""

from __future__ import annotations

import numpy as np

from outputs import BY_NAME, FIELDS, reference_sequences, shape

# Small enough to write out by hand and to read in a failure, and ragged enough
# that a row, a histogram and a scalar are all of different widths.
N_REFS = 4
CAP = 6

# A count past what a float32 holds exactly, so a value that comes back whole
# was not narrowed to the type the rates use.
LARGE_COUNT = np.uint64(2) ** 40 + 1

# Text for a file that is not an output. A run that fails must not modify it.
NOTES = "months of irreplaceable notes\n"


def not_hdf5(tmp_path):
    """Writes a file that is not an output, for a run to fail on."""
    notes = tmp_path / "notes.txt"
    notes.write_text(NOTES)

    return notes


def random_values(field, n_refs=N_REFS, cap=CAP, seed=0):
    """Builds values for one field that differ from column to column and from
    row to row, so that a rule applied along the wrong axis does not match by
    accident."""
    rng = np.random.default_rng(seed)
    wanted = shape(BY_NAME[field], n_refs, cap)

    if BY_NAME[field].sequence:
        return reference_sequences(n_refs, cap)

    if BY_NAME[field].dtype == "u8":
        return rng.integers(0, 1000, size=wanted, dtype=np.uint64)

    return rng.random(size=wanted).astype(np.float32)


def random_fields(seed, n_refs=N_REFS, cap=CAP) -> dict:
    """Builds values for every field at once, so that a test of one field runs
    on a file whose other fields are not zero."""
    return {field.name: random_values(field.name, n_refs, cap, seed=seed + i)
            for i, field in enumerate(FIELDS)}


def missing_in_each_input(rows: int, columns: int):
    """Builds a pair of arrays that are NaN in the first input only, in the
    second input only, and in both."""
    known = np.float32(0.5)

    left = np.full((rows, columns), known, dtype=np.float32)
    right = np.full((rows, columns), known, dtype=np.float32)

    left[1, :] = np.nan
    right[2, :] = np.nan
    left[3, :] = right[3, :] = np.nan

    return left, right
