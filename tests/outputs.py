"""Building output files by hand, and the arithmetic over two of them.

cmuts-sub takes two outputs and writes a third. Its result depends only on the
values in those files and not on how they were produced, so nothing here runs
cmuts to obtain them: the inputs are written directly, with values chosen to
exercise each rule. Tests built this way cover the subtraction alone, and remain
valid if the reactivity calculation changes.

The layout is the one thing the two programs share. It is written out here
rather than read from a cmuts run, so that a field added to the output fails a
single named test rather than the whole suite; see test_subtraction.py, which
checks this description against a real run.
"""

from __future__ import annotations

from dataclasses import dataclass

import h5py
import numpy as np

# One value per reference base, one per read length, one per reference.
PER_BASE = "per base"
PER_LENGTH = "per length"
SCALAR = "scalar"


@dataclass(frozen=True)
class Field:
    name: str
    kind: str
    dtype: str
    fill: float


# The datasets an output holds, as output.h declares them. A rate is float and
# missing where it was not measured; a count is a whole unsigned and zero where
# nothing was counted.
FIELDS = (
    Field("coverage", PER_BASE, "f4", 0.0),
    Field("reactivity", PER_BASE, "f4", np.nan),
    Field("error", PER_BASE, "f4", np.nan),
    Field("reads/lengths", PER_LENGTH, "u8", 0),
    Field("reads/counted", SCALAR, "u8", 0),
    Field("reads/rejected", SCALAR, "u8", 0),
)

# Belongs to the run and not to any reference, so it has no row.
RUN_TOTAL = "reads/unmapped"

READS_GROUP = "reads"

BY_NAME = {field.name: field for field in FIELDS}


# ---------------------------------------------------------------------------
# The rules
# ---------------------------------------------------------------------------


def _add(treated, untreated):
    return treated + untreated


def _subtract(treated, untreated):
    return treated - untreated


def _quadrature(treated, untreated):
    return np.sqrt(treated * treated + untreated * untreated)


# How each dataset of the output is formed from the two inputs. Stated here as
# the contract, and not derived from the program's own description of it.
RULES = {
    "coverage": _add,
    "reactivity": _subtract,
    "error": _quadrature,
    "reads/lengths": _add,
    "reads/counted": _add,
    "reads/rejected": _add,
    RUN_TOTAL: _add,
}


def combined(name, treated, untreated):
    """What the output should hold, given what the two inputs hold.

    The rates are stored as float32 and combined as double, so the expectation
    is formed the same way round: widened, combined, narrowed. Comparisons
    against it are therefore exact, and a combination carried out in float32
    would not match.
    """
    rule = RULES[name]
    treated, untreated = np.asarray(treated), np.asarray(untreated)

    if np.issubdtype(treated.dtype, np.floating):
        return rule(treated.astype(np.float64),
                    untreated.astype(np.float64)).astype(np.float32)

    return rule(treated, untreated)


# ---------------------------------------------------------------------------
# Writing one
# ---------------------------------------------------------------------------


def width(kind: str, cap: int) -> int:
    """Columns a field of this kind occupies, for a file whose longest
    reference is cap bases.

    A read length is not a position in a reference: a read carrying insertions
    or soft-clipped ends is longer than what it aligns to, so the histogram
    reaches to twice the longest reference. It starts at length 1, a read
    storing no sequence never reaching the tally, so bin i holds length i + 1.
    """
    return cap if kind == PER_BASE else 2 * cap


def shape(field: Field, n_refs: int, cap: int) -> tuple:
    return (n_refs,) if field.kind == SCALAR else (n_refs, width(field.kind, cap))


def _values(field: Field, n_refs: int, cap: int, given) -> np.ndarray:
    """A field's dataset, as the caller specified it: an array of the right
    shape, or a single value filling the whole of it."""
    wanted = shape(field, n_refs, cap)

    if given is None:
        given = field.fill

    if np.isscalar(given):
        return np.full(wanted, given, dtype=field.dtype)

    return np.asarray(given, dtype=field.dtype).reshape(wanted)


def _storage(storage: str, wanted: tuple) -> dict:
    """How the datasets are stored, which a reader's result must not depend on.
    cmuts writes chunked, shuffled and deflated; a plainly written file holds
    the same values and is the other case worth reading.

    A dataset with no rows cannot be chunked, so an empty one is written plainly
    whichever was asked for.
    """
    if storage == "plain" or 0 in wanted:
        return {}

    return dict(chunks=wanted, shuffle=True, compression="gzip", compression_opts=3)


def write_output(path, *, n_refs=4, cap=6, values=None, unmapped=0, storage="plain"):
    """Writes a file in the output layout. A field the caller does not specify
    is filled with the value marking it as never written."""
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

        out.create_dataset(RUN_TOTAL, data=np.uint64(unmapped))

    return path


# ---------------------------------------------------------------------------
# Reading one, and damaging one
# ---------------------------------------------------------------------------


def field_of(path, name) -> np.ndarray:
    """One dataset, read whole. [()] rather than [:], a run total being a
    scalar that cannot be sliced."""
    with h5py.File(path, "r") as handle:
        return handle[name][()]


def datasets_of(path) -> dict:
    """Every dataset of a file, by path, with its shape and dtype. The counts
    about reads sit in a group of their own, so a name is a path and not a
    key."""
    found = {}

    with h5py.File(path, "r") as handle:
        handle.visititems(
            lambda name, obj: found.update({name: (obj.shape, obj.dtype)})
            if isinstance(obj, h5py.Dataset) else None
        )

    return found


def without(path, name):
    """The same file with one object removed, for the inputs that should be
    refused."""
    with h5py.File(path, "r+") as handle:
        del handle[name]

    return path


def rewidened(path, name, columns):
    """The same file with one dataset given a width the rest disagree with."""
    with h5py.File(path, "r+") as handle:
        rows = handle[name].shape[0]
        dtype = handle[name].dtype
        del handle[name]
        handle.create_dataset(name, shape=(rows, columns), dtype=dtype)

    return path
