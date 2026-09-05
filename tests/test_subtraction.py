"""Taking an untreated background off a treated run.

The result depends only on the values in the input files, so the
inputs are written by hand and not counted from an alignment. inputs.py builds
them and outputs.py describes the layout the programs share. The contracts
this program shares with the other readers of outputs are in test_io.py.
"""

from __future__ import annotations

import numpy as np
import pytest

from inputs import CAP, N_REFS, missing_in_each_input, random_fields, random_values
from outputs import ALL_FIELDS, ERROR_FIELDS, RATE_FIELDS, field_of
from programs import run_subtract

rates = pytest.mark.parametrize("rate", RATE_FIELDS)
errors = pytest.mark.parametrize("error", ERROR_FIELDS)


@pytest.fixture
def subtract(tmp_path):
    """Returns a function that runs the subtraction into a path of its own."""
    def run(treated, untreated, **options):
        return run_subtract(treated, untreated, tmp_path / "difference.h5", **options)

    return run


@rates
def test_a_background_above_the_signal_leaves_a_negative_rate(build, subtract, rate):
    treated = build({rate: 0.25})
    untreated = build({rate: 0.75})

    difference = field_of(subtract(treated, untreated), rate)

    assert np.all(difference == np.float32(-0.5))


@rates
def test_clipping_holds_the_difference_at_zero(build, subtract, rate):
    treated = build({rate: 0.25})
    untreated = build({rate: 0.75})

    output = subtract(treated, untreated, clip=True)

    assert np.all(field_of(output, rate) == 0)


@rates
def test_clipping_leaves_a_difference_above_zero_alone(build, tmp_path, rate):
    treated, untreated = build(random_fields(seed=13)), build(random_fields(seed=14))

    plain = run_subtract(treated, untreated, tmp_path / "plain.h5")
    clipped = run_subtract(treated, untreated, tmp_path / "clipped.h5", clip=True)

    unclipped = field_of(plain, rate)

    assert (unclipped < 0).any()
    assert np.array_equal(field_of(clipped, rate), np.maximum(unclipped, 0))


@rates
def test_clipping_does_not_raise_a_missing_value_to_zero(build, subtract, rate):
    left, right = missing_in_each_input(N_REFS, CAP)

    output = subtract(build({rate: left}), build({rate: right}), clip=True)
    result = field_of(output, rate)

    assert np.array_equal(np.isnan(result), np.isnan(left) | np.isnan(right))


def test_clipping_reaches_no_field_but_the_rates(build, tmp_path):
    treated = build(random_fields(seed=15), unmapped=11)
    untreated = build(random_fields(seed=16), unmapped=4)

    plain = run_subtract(treated, untreated, tmp_path / "plain.h5")
    clipped = run_subtract(treated, untreated, tmp_path / "clipped.h5", clip=True)

    for name in set(ALL_FIELDS) - set(RATE_FIELDS):
        assert np.array_equal(field_of(clipped, name), field_of(plain, name)), name


@errors
def test_an_error_is_never_reduced_by_subtracting(build, subtract, error):
    treated = build({error: random_values(error, seed=7)})
    untreated = build({error: random_values(error, seed=8)})

    output = subtract(treated, untreated)

    assert np.all(field_of(output, error) >= field_of(treated, error))
    assert np.all(field_of(output, error) >= field_of(untreated, error))


@rates
def test_a_file_against_itself_leaves_a_rate_of_zero(build, subtract, rate):
    values = random_fields(seed=3)
    values[rate][2, 3] = np.nan

    treated = build(values)
    result = field_of(subtract(treated, treated), rate)
    known = ~np.isnan(field_of(treated, rate))

    assert known.any()
    assert np.all(result[known] == 0)
    assert np.isnan(result[~known]).all()


def test_swapping_the_two_inputs_negates_only_the_rates(build, subtract, tmp_path):
    treated, untreated = build(random_fields(seed=4)), build(random_fields(seed=40))

    forward = subtract(treated, untreated)
    backward = run_subtract(untreated, treated, tmp_path / "backward.h5")

    for rate in RATE_FIELDS:
        assert np.array_equal(field_of(backward, rate),
                              -field_of(forward, rate)), rate

    for name in set(ALL_FIELDS) - set(RATE_FIELDS):
        assert np.array_equal(field_of(backward, name), field_of(forward, name)), name


def test_a_background_of_zeros_leaves_the_treated_run_unchanged(build, subtract):
    treated = build(random_fields(seed=5))
    zeros = build(dict.fromkeys(RATE_FIELDS + ERROR_FIELDS, 0.0))

    output = subtract(treated, zeros)

    for name in ALL_FIELDS:
        assert np.array_equal(field_of(output, name), field_of(treated, name)), name
