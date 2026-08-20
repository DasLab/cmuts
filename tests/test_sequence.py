"""The sequence dataset, which holds the reference itself.
"""

import h5py
import numpy as np
import pytest

from oracle import sequences
from outputs import REACTIVITY, SEQUENCE, delete_field, field_of, layout_of
from programs import (
    run_cmuts,
    run_divide,
    run_normalize,
    run_subtract,
    try_divide,
    try_subtract,
)

# What each base is written as. U and T are one base, and anything else is
# written as the marker a column outside a reference carries.
CODE = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}
OUTSIDE = -1


def encoded(bases: str, width: int) -> np.ndarray:
    """Returns the row a reference of these bases should occupy."""
    row = np.full(width, OUTSIDE, dtype=np.int8)

    for i, base in enumerate(bases.upper()):
        row[i] = CODE.get(base, OUTSIDE)

    return row


def alter_a_base(path):
    """Changes one base of one reference, so the file describes a library no
    other input describes."""
    with h5py.File(path, "r+") as handle:
        handle[SEQUENCE][0, 0] = (handle[SEQUENCE][0, 0] + 1) % 4

    return path


# ---------------------------------------------------------------------------
# What cmuts hmm writes
# ---------------------------------------------------------------------------


def test_every_reference_holds_its_own_bases(data, tmp_path, falsifiable):
    output = tmp_path / "counted.h5"
    run_cmuts(data, output)

    written = field_of(output, SEQUENCE)
    fasta = sequences(data.fasta)

    falsifiable(len(fasta) > 0)

    assert written.dtype == np.int8
    assert written.shape[0] == len(fasta)

    for row, (name, bases) in enumerate(fasta.items()):
        assert np.array_equal(written[row], encoded(bases, written.shape[1])), name


# ---------------------------------------------------------------------------
# What the programs reading an output do with it
# ---------------------------------------------------------------------------


def test_subtraction_carries_the_sequence_through(build, tmp_path):
    treated, untreated = build(), build()
    output = run_subtract(treated, untreated, tmp_path / "difference.h5")

    assert np.array_equal(field_of(output, SEQUENCE), field_of(treated, SEQUENCE))


def test_division_carries_the_sequence_through(build, tmp_path):
    rates, control = build(), build()
    output = run_divide(rates, control, tmp_path / "ratio.h5")

    assert np.array_equal(field_of(output, SEQUENCE), field_of(rates, SEQUENCE))


def test_normalization_carries_the_sequence_through(build, tmp_path):
    rates = build()
    run_normalize([rates], [tmp_path / "norm.h5"])

    assert np.array_equal(field_of(tmp_path / "norm.h5", SEQUENCE),
                          field_of(rates, SEQUENCE))


@pytest.mark.parametrize("attempt", [try_subtract, try_divide])
def test_inputs_whose_sequences_differ_are_refused(build, tmp_path, attempt):
    failed = attempt(build(), alter_a_base(build()), tmp_path / "out.h5")

    assert failed.returncode != 0
    assert SEQUENCE in failed.stderr


# ---------------------------------------------------------------------------
# An output written before the dataset existed
# ---------------------------------------------------------------------------


def test_subtraction_omits_the_sequence_when_an_input_has_none(build, tmp_path):
    older = delete_field(build(), SEQUENCE)
    output = run_subtract(build(), older, tmp_path / "difference.h5")

    assert SEQUENCE not in layout_of(output)
    assert REACTIVITY in layout_of(output)


def test_normalization_omits_the_sequence_when_the_input_has_none(build, tmp_path):
    older = delete_field(build(), SEQUENCE)
    run_normalize([older], [tmp_path / "norm.h5"])

    assert SEQUENCE not in layout_of(tmp_path / "norm.h5")
    assert REACTIVITY in layout_of(tmp_path / "norm.h5")
