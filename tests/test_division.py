"""Dividing the rates by a denatured control.

The result depends only on the values in the input files, so the
inputs are written by hand and not counted from an alignment. inputs.py builds
them and outputs.py describes the layout the programs share. The contracts
this program shares with the other readers of outputs are in test_io.py.
"""

from __future__ import annotations

import numpy as np
import pytest

from inputs import CAP, N_REFS, random_fields, random_values
from outputs import COUNTED, ERROR_FIELDS, RATE_FIELDS, UNMAPPED, field_of
from programs import run_divide, run_subtract

# The error rounds at every step of the division and the root, so it is compared
# to a relative tolerance. Every other rule rounds once and is exact.
TOLERANCE = 1e-6

# One id per channel, taken from the group half of each rate's path.
CHANNELS = tuple(name.split("/")[0] for name in RATE_FIELDS)

channels = pytest.mark.parametrize("rate, error", tuple(zip(RATE_FIELDS, ERROR_FIELDS)),
                                   ids=CHANNELS)
rates = pytest.mark.parametrize("rate", RATE_FIELDS, ids=CHANNELS)


@pytest.fixture
def divide(tmp_path):
    """Returns a function that runs the division into a path of its own."""
    def run(rates, control, **options):
        return run_divide(rates, control, tmp_path / "normalized.h5", **options)

    return run


# ---------------------------------------------------------------------------
# The control
# ---------------------------------------------------------------------------


@rates
def test_a_control_divides_the_rates(build, divide, rate):
    treated = build({rate: 0.5})
    control = build({rate: 0.25})

    output = divide(treated, control)

    assert np.all(field_of(output, rate) == np.float32(2.0))


def test_a_control_of_ones_leaves_the_rates_alone(build, divide):
    treated = build(random_fields(seed=34))
    ones = build(dict.fromkeys(RATE_FIELDS, 1.0) | dict.fromkeys(ERROR_FIELDS, 0.0))

    output = divide(treated, ones)

    for name in RATE_FIELDS + ERROR_FIELDS:
        assert np.array_equal(field_of(output, name), field_of(treated, name)), name


@channels
@pytest.mark.parametrize("value", [0.0, -0.25])
def test_a_control_not_above_zero_leaves_no_rate(build, divide, rate, error, value):
    output = divide(build({rate: 0.5}), build({rate: value}))

    assert np.isnan(field_of(output, rate)).all()
    assert np.isnan(field_of(output, error)).all()


@rates
def test_a_control_of_nan_leaves_no_rate(build, divide, rate):
    control = np.full((N_REFS, CAP), np.float32(0.5), dtype=np.float32)
    control[2, :] = np.nan

    output = divide(build({rate: 0.75}), build({rate: control}))
    result = field_of(output, rate)

    assert np.isnan(result[2]).all()
    assert not np.isnan(result[0]).any()


@channels
def test_an_uncertain_control_widens_the_error(build, tmp_path, rate, error):
    treated = build(random_fields(seed=36))
    exact = build({rate: 1.0, error: 0.0})
    uncertain = build({rate: 1.0, error: random_values(error, seed=38)})

    against_exact = run_divide(treated, exact, tmp_path / "exact.h5")
    against_uncertain = run_divide(treated, uncertain, tmp_path / "uncertain.h5")

    assert np.all(field_of(against_uncertain, error) >= field_of(against_exact, error))


def test_a_control_is_counted_in_the_totals(build, divide):
    treated = build({COUNTED: 3}, unmapped=11)
    control = build({COUNTED: 5}, unmapped=13)

    output = divide(treated, control)

    assert np.all(field_of(output, COUNTED) == 8)
    assert field_of(output, UNMAPPED) == 24


@rates
def test_a_file_against_itself_leaves_a_rate_of_one(build, divide, rate):
    values = random_fields(seed=3)
    values[rate] = np.abs(values[rate]) + np.float32(0.25)
    values[rate][2, 3] = np.nan

    treated = build(values)
    result = field_of(divide(treated, treated), rate)
    known = ~np.isnan(field_of(treated, rate))

    assert known.any()
    assert np.all(result[known] == np.float32(1.0))
    assert np.isnan(result[~known]).all()


# ---------------------------------------------------------------------------
# A subtraction followed by a division
# ---------------------------------------------------------------------------


def _normalized(tmp_path, treated, untreated, control, **options):
    """Subtracts the background and divides the difference by the control."""
    difference = run_subtract(treated, untreated, tmp_path / "difference.h5", **options)

    return run_divide(difference, control, tmp_path / "normalized.h5")


@rates
def test_the_two_programs_give_the_rate_over_the_control(build, tmp_path, rate):
    treated = build(random_fields(seed=41))
    untreated = build(random_fields(seed=42))
    control = build(random_fields(seed=43))

    output = _normalized(tmp_path, treated, untreated, control)

    treated_rate = field_of(treated, rate)
    untreated_rate = field_of(untreated, rate)
    control_rate = field_of(control, rate)

    with np.errstate(divide="ignore", invalid="ignore"):
        normalized = np.where(control_rate > 0,
                              (treated_rate - untreated_rate) / control_rate,
                              np.float32(np.nan))

    assert np.array_equal(field_of(output, rate), normalized, equal_nan=True)


@channels
def test_the_two_programs_add_the_errors_in_quadrature(build, tmp_path, rate, error):
    treated = build(random_fields(seed=44))
    untreated = build(random_fields(seed=45))
    control = build(random_fields(seed=46))

    output = _normalized(tmp_path, treated, untreated, control)

    treated_error = field_of(treated, error)
    untreated_error = field_of(untreated, error)
    control_error = field_of(control, error)
    control_rate = field_of(control, rate)
    normalized_rate = field_of(output, rate)

    variance = (treated_error ** 2 + untreated_error ** 2
                + (normalized_rate * control_error) ** 2)

    with np.errstate(divide="ignore", invalid="ignore"):
        normalized = np.where(control_rate > 0, np.sqrt(variance) / control_rate,
                              np.float32(np.nan))

    assert np.allclose(field_of(output, error), normalized, rtol=TOLERANCE,
                       equal_nan=True)


def test_every_count_is_summed_over_all_three_inputs(build, tmp_path):
    treated = build({COUNTED: 3}, unmapped=11)
    untreated = build({COUNTED: 5}, unmapped=13)
    control = build({COUNTED: 7}, unmapped=17)

    output = _normalized(tmp_path, treated, untreated, control)

    assert np.all(field_of(output, COUNTED) == 15)
    assert field_of(output, UNMAPPED) == 41


@rates
def test_clipping_the_difference_holds_the_normalized_rate_at_zero(build, tmp_path,
                                                                   rate):
    treated = build({rate: 0.25})
    untreated = build({rate: 0.75})
    control = build({rate: 0.5})

    output = _normalized(tmp_path, treated, untreated, control, clip=True)

    assert np.all(field_of(output, rate) == 0)


@channels
def test_a_clipped_rate_carries_no_uncertainty_from_the_control(build, tmp_path,
                                                                rate, error):
    """The error of a ratio scales the control's error by the rate, so at a rate
    of zero the control contributes no term and the error is the quadrature of
    the two inputs divided by the control's rate."""
    treated = build({rate: 0.25, error: 0.3})
    untreated = build({rate: 0.75, error: 0.4})
    control = build({rate: 0.5, error: 0.1})

    output = _normalized(tmp_path, treated, untreated, control, clip=True)

    quadrature = np.sqrt(np.float32(0.3) ** 2 + np.float32(0.4) ** 2)

    assert np.allclose(field_of(output, error), quadrature / np.float32(0.5),
                       rtol=TOLERANCE)
