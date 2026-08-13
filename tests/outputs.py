"""Building output files by hand, and the arithmetic over two of them.

cmuts-sub takes two outputs and writes a third. Its result depends only on the
values in those files and not on how they were produced, so nothing here runs
cmuts to obtain them: the inputs are written directly, with values chosen to
exercise each rule. Tests built this way cover only the subtraction, and remain
valid if the reactivity calculation changes.

The layout is the one thing the two programs share. It is written out here
rather than read from a cmuts run, so that the description is a contract and
not a copy of whatever cmuts currently produces.
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
    rate: bool = False


# The datasets an output holds, as output.h declares them. A rate is float and
# missing where it was not measured; a count is a whole unsigned and zero where
# nothing was counted. A rate is also the one kind two files hold one value of
# rather than a total, so nothing adds it up.
FIELDS = (
    Field("coverage", PER_BASE, "f4", 0.0),
    Field("reactivity", PER_BASE, "f4", np.nan, rate=True),
    Field("error", PER_BASE, "f4", np.nan, rate=True),
    Field("reads/lengths", PER_LENGTH, "u8", 0),
    Field("reads/counted", SCALAR, "u8", 0),
    Field("reads/rejected", SCALAR, "u8", 0),
)

# Belongs to the run and not to any reference, so it has no row.
RUN_TOTAL = "reads/unmapped"

READS_GROUP = "reads"

BY_NAME = {field.name: field for field in FIELDS}

# The groups the tests sort the datasets into, each read off the declaration
# above rather than listed again beside it.
RATES = tuple(field.name for field in FIELDS if field.rate)
BY_LENGTH = tuple(field.name for field in FIELDS if field.kind == PER_LENGTH)
FLOATING = tuple(field.name for field in FIELDS if field.dtype.startswith("f"))


# ---------------------------------------------------------------------------
# The rules
# ---------------------------------------------------------------------------


# Every rule takes the values of the field it forms and the reactivities of the
# same inputs, one array apiece. Most ignore the second: only the error of a
# ratio reads a field other than its own.


def _add(values, rates):
    total = values[0]

    for value in values[1:]:
        total = total + value

    return total


def _subtract(values, rates):
    return values[0] - values[1]


def _quadrature(values, rates):
    return np.sqrt(values[0] * values[0] + values[1] * values[1])


def _over_control(values, denatured):
    """A rate over the denatured control's, which is undefined where the control
    measured nothing to divide by."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(denatured > 0, values / denatured, np.float32(np.nan))


def _ratio(values, rates):
    return _over_control(values[0] - values[1], values[2])


def _ratio_error(values, rates):
    ratio = _over_control(rates[0] - rates[1], rates[2])
    spread = values[0] ** 2 + values[1] ** 2 + (ratio * values[2]) ** 2

    return _over_control(np.sqrt(spread), rates[2])


# How each dataset of the output is formed from the inputs, without a denatured
# control and with one. Stated here as the contract, and not derived from the
# program's own description of it.
RULES = {
    "coverage": (_add, _add),
    "reactivity": (_subtract, _ratio),
    "error": (_quadrature, _ratio_error),
    "reads/lengths": (_add, _add),
    "reads/counted": (_add, _add),
    "reads/rejected": (_add, _add),
    RUN_TOTAL: (_add, _add),
}

CONTROLLED = 1
UNCONTROLLED = 0


def combined(name, *inputs):
    """What the output should hold, given the files that went into it.

    Takes the files rather than one field's values from each, since the error
    of a ratio reads the reactivities as well as the errors.

    Each rule is applied in the type the field is stored as, which is what
    cmuts-sub does: rates in float32, counts as whole unsigneds. Comparisons
    against a sum or a difference are therefore exact -- widening to double
    first would agree, one rounding standing in for the other -- but a ratio
    and its error round often enough that a caller should allow a tolerance.
    """
    rule = RULES[name][CONTROLLED if len(inputs) > 2 else UNCONTROLLED]
    values = [np.asarray(field_of(path, name)) for path in inputs]
    rates = [np.asarray(field_of(path, "reactivity")) for path in inputs]

    return rule(values, rates)


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

    A dataset with no rows cannot be chunked, so this writes an empty one
    plainly under either setting.
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


def layout_of(path) -> dict:
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
