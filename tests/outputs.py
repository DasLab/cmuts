"""Tools for reading and writing HDF5 files in the format cmuts uses."""

from __future__ import annotations

from dataclasses import dataclass

import h5py
import numpy as np

COVERAGE = "coverage"
REACTIVITY = "reactivity"
ERROR = "error"
LENGTHS = "reads/lengths"
COUNTED = "reads/counted"
REJECTED = "reads/rejected"
UNMAPPED = "reads/unmapped"
NORM = "norm"

READS_GROUP = "reads"

# One value per reference base, one per read length, one per reference.
PER_BASE = "per base"
PER_LENGTH = "per length"
SCALAR = "scalar"


@dataclass(frozen=True)
class Field:
    """A dataset of the output, as declared in output.h."""

    name: str
    kind: str
    dtype: str
    fill: float          # what a row nobody wrote reads as
    pad: float           # what a column past the end of a reference reads as
    rate: bool = False


# A rate is a float and NaN where it was not measured; a count is a whole
# unsigned and zero where nothing was counted.
FIELDS = (
    Field(COVERAGE, PER_BASE, "f4", 0.0, np.nan),
    Field(REACTIVITY, PER_BASE, "f4", np.nan, np.nan, rate=True),
    Field(ERROR, PER_BASE, "f4", np.nan, np.nan, rate=True),
    Field(LENGTHS, PER_LENGTH, "u8", 0, 0),
    Field(COUNTED, SCALAR, "u8", 0, 0),
    Field(REJECTED, SCALAR, "u8", 0, 0),
)

BY_NAME = {field.name: field for field in FIELDS}


def _field_names(matches) -> tuple:
    """Returns the names of the fields matching a predicate."""
    return tuple(field.name for field in FIELDS if matches(field))


# The groups a test reasons about. Each is derived from the declaration above,
# so a new field joins them without being listed a second time.
PER_BASE_FIELDS = _field_names(lambda field: field.kind == PER_BASE)
PER_LENGTH_FIELDS = _field_names(lambda field: field.kind == PER_LENGTH)
PER_REFERENCE_FIELDS = _field_names(lambda field: field.kind == SCALAR)
RATE_FIELDS = _field_names(lambda field: field.rate)
FLOAT_FIELDS = _field_names(lambda field: field.dtype.startswith("f"))

# The fields holding counts, in which a zero is a measured value and not a
# missing one.
COUNT_FIELDS = _field_names(lambda field: not field.rate)
COLUMN_COUNT_FIELDS = _field_names(lambda field: not field.rate and field.kind != SCALAR)

# The fields in which a NaN marks a column past the end of a reference. A row
# indexed by read length runs to its full width, and a rate is NaN at any
# position failing --min-depth, so neither kind marks padding.
PADDED_FIELDS = _field_names(lambda field: not field.rate and field.kind == PER_BASE)

# Every field with a row per reference, and every field of the file. The
# unmapped total belongs to the run and not to any reference, so it has no row
# and joins none of the groups above.
ROW_FIELDS = tuple(field.name for field in FIELDS)
ALL_FIELDS = ROW_FIELDS + (UNMAPPED,)


def width(kind: str, cap: int) -> int:
    """Returns the columns a field of this kind occupies, for a file whose
    longest reference is cap bases.

    Insertions and soft-clipped ends make a read longer than what it aligns to,
    so the length histogram reaches to twice the longest reference. It starts
    at length 1, so bin i holds length i + 1.
    """
    return cap if kind == PER_BASE else 2 * cap


def shape(field: Field, n_refs: int, cap: int) -> tuple:
    """Returns the shape a field takes in a file of this size."""
    return (n_refs,) if field.kind == SCALAR else (n_refs, width(field.kind, cap))


def extent(name: str, reference_length: int, columns: int) -> int:
    """Returns how far along a row the field holds data; the rest is padding."""
    return columns if name in PER_LENGTH_FIELDS else reference_length


# ---------------------------------------------------------------------------
# Writing an output
# ---------------------------------------------------------------------------


def _values(field: Field, n_refs: int, cap: int, given) -> np.ndarray:
    """Builds a field's dataset from what the caller specified: an array of the
    right shape, or a single value to fill it with."""
    wanted = shape(field, n_refs, cap)

    if given is None:
        given = field.fill

    if np.isscalar(given):
        return np.full(wanted, given, dtype=field.dtype)

    return np.asarray(given, dtype=field.dtype).reshape(wanted)


def _storage(storage: str, wanted: tuple) -> dict:
    """Returns the creation settings a dataset is written under. cmuts hmm
    writes chunked, shuffled and deflated, and a plainly written file holds the
    same values.

    A dataset with no rows cannot be chunked, so an empty one is written plainly
    under either setting.
    """
    if storage == "plain" or 0 in wanted:
        return {}

    return dict(chunks=wanted, shuffle=True, compression="gzip", compression_opts=3)


def write_output(path, *, n_refs=4, cap=6, values=None, unmapped=0, storage="plain"):
    """Writes a file in the output layout. A field the caller does not specify
    is filled with its fill value: NaN for a rate, zero for a count."""
    values = values or {}
    unknown = set(values) - set(BY_NAME)
    if unknown:
        raise KeyError(f"no such field: {', '.join(sorted(unknown))}")

    with h5py.File(path, "w") as out:
        out.create_group(READS_GROUP)

        for field in FIELDS:
            data = _values(field, n_refs, cap, values.get(field.name))
            out.create_dataset(field.name, data=data,
                               **_storage(storage, data.shape))

        out.create_dataset(UNMAPPED, data=np.uint64(unmapped))

    return path


# ---------------------------------------------------------------------------
# Reading an output
# ---------------------------------------------------------------------------


def arrays_of(handle) -> dict:
    """Returns every dataset of an open file, by path. The counts about reads
    sit in a group of their own, so a name is a path and not a key."""
    found = {}

    handle.visititems(
        lambda name, obj: found.update({name: obj})
        if isinstance(obj, h5py.Dataset) else None
    )

    return found


def field_of(path, name) -> np.ndarray:
    """Reads a single dataset in full. A run total is a scalar, so [()] is used
    to index it."""
    with h5py.File(path, "r") as handle:
        return handle[name][()]


def fields_of(path, names) -> dict:
    """Reads the named datasets, keyed by name."""
    with h5py.File(path, "r") as handle:
        return {name: handle[name][()] for name in names}


def attributes_of(path) -> dict:
    """Reads the attributes of the root group, as text. HDF5 hands a string back
    as bytes, whatever it was written from."""
    with h5py.File(path, "r") as handle:
        return {name: value.decode() if isinstance(value, bytes) else value
                for name, value in handle.attrs.items()}


def layout_of(path) -> dict:
    """Returns the shape and type of every dataset, by path."""
    with h5py.File(path, "r") as handle:
        return {name: (array.shape, array.dtype)
                for name, array in arrays_of(handle).items()}


@dataclass(frozen=True)
class Summary:
    """The read counts a run wrote."""

    kept: int
    rejected: int
    rows: int
    unmapped: int


def read_summary(path) -> Summary:
    """Reads the counts a run made over the reads it was given."""
    with h5py.File(path, "r") as output:
        counted = output[COUNTED][:]
        rejected = output[REJECTED][:]

        return Summary(
            kept=int(np.nansum(counted, dtype=np.float64)),
            rejected=int(np.nansum(rejected, dtype=np.float64)),
            # The counts are zero-filled, so a row with a nonzero total is one
            # that at least one read aligned to, kept or rejected.
            rows=int(((counted + rejected) > 0).sum()),
            unmapped=int(output[UNMAPPED][()]),
        )


def arrays_agree(a, b) -> bool:
    """Returns whether two arrays hold the same values. Two floating ones agree where
    both hold a NaN in the same position."""
    return np.array_equal(a[()], b[()], equal_nan=np.issubdtype(a.dtype, np.floating))


def outputs_agree(first, second) -> bool:
    """Returns whether two files hold the same thing.

    Compares every dataset the files hold, so a new field is checked as soon as
    it is added.
    """
    with h5py.File(first, "r") as a, h5py.File(second, "r") as b:
        left, right = arrays_of(a), arrays_of(b)

        return set(left) == set(right) and all(
            arrays_agree(left[name], right[name]) for name in left
        )


# ---------------------------------------------------------------------------
# Producing malformed outputs
# ---------------------------------------------------------------------------


def delete_field(path, name):
    """Removes a dataset from a file."""
    with h5py.File(path, "r+") as handle:
        del handle[name]

    return path


def set_field_width(path, name, columns):
    """Sets the column count of a dataset, leaving it inconsistent with the
    rest of the file."""
    with h5py.File(path, "r+") as handle:
        rows = handle[name].shape[0]
        dtype = handle[name].dtype
        del handle[name]
        handle.create_dataset(name, shape=(rows, columns), dtype=dtype)

    return path
