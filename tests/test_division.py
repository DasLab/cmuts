"""Dividing reactivity rates by a denatured control.

The result depends only on the values in the input files, so the
inputs are written by hand and not counted from an alignment. inputs.py builds
them and outputs.py describes the layout the programs share. The contracts
this program shares with the other readers of outputs are in test_io.py.
"""

from __future__ import annotations

import numpy as np
import pytest

from inputs import CAP, N_REFS, random_fields, random_values
from outputs import COUNTED, ERROR, REACTIVITY, UNMAPPED, field_of
from programs import run_divide, run_subtract

# The error rounds at every step of the division and the root, so it is compared
# to a relative tolerance. Every other rule rounds once and is exact.
TOLERANCE = 1e-6


@pytest.fixture
def divide(tmp_path):
    """Returns a function that runs the division into a path of its own."""
    def run(rates, control, **options):
        return run_divide(rates, control, tmp_path / "normalized.h5", **options)

    return run


# ---------------------------------------------------------------------------
# The control
# ---------------------------------------------------------------------------


def test_a_control_divides_the_rates(build, divide):
    rates = build({REACTIVITY: 0.5})
    control = build({REACTIVITY: 0.25})

    output = divide(rates, control)

    assert np.all(field_of(output, REACTIVITY) == np.float32(2.0))


def test_a_control_of_ones_leaves_the_rates_alone(build, divide):
    rates = build(random_fields(seed=34))
    ones = build({REACTIVITY: 1.0, ERROR: 0.0})

    output = divide(rates, ones)

    assert np.array_equal(field_of(output, REACTIVITY), field_of(rates, REACTIVITY))
    assert np.array_equal(field_of(output, ERROR), field_of(rates, ERROR))


@pytest.mark.parametrize("rate", [0.0, -0.25])
def test_a_control_not_above_zero_leaves_no_reactivity(build, divide, rate):
    output = divide(build({REACTIVITY: 0.5}), build({REACTIVITY: rate}))

    assert np.isnan(field_of(output, REACTIVITY)).all()
    assert np.isnan(field_of(output, ERROR)).all()


def test_a_control_of_nan_leaves_no_reactivity(build, divide):
    control = np.full((N_REFS, CAP), np.float32(0.5), dtype=np.float32)
    control[2, :] = np.nan

    output = divide(build({REACTIVITY: 0.75}), build({REACTIVITY: control}))
    result = field_of(output, REACTIVITY)

    assert np.isnan(result[2]).all()
    assert not np.isnan(result[0]).any()


def test_an_uncertain_control_widens_the_error(build, tmp_path):
    rates = build(random_fields(seed=36))
    exact = build({REACTIVITY: 1.0, ERROR: 0.0})
    uncertain = build({REACTIVITY: 1.0, ERROR: random_values(ERROR, seed=38)})

    against_exact = run_divide(rates, exact, tmp_path / "exact.h5")
    against_uncertain = run_divide(rates, uncertain, tmp_path / "uncertain.h5")

    assert np.all(field_of(against_uncertain, ERROR) >= field_of(against_exact, ERROR))


def test_a_control_is_counted_in_the_totals(build, divide):
    rates = build({COUNTED: 3}, unmapped=11)
    control = build({COUNTED: 5}, unmapped=13)

    output = divide(rates, control)

    assert np.all(field_of(output, COUNTED) == 8)
    assert field_of(output, UNMAPPED) == 24


def test_a_file_against_itself_leaves_a_reactivity_of_one(build, divide):
    values = random_fields(seed=3)
    values[REACTIVITY] = np.abs(values[REACTIVITY]) + np.float32(0.25)
    values[REACTIVITY][2, 3] = np.nan

    rates = build(values)
    result = field_of(divide(rates, rates), REACTIVITY)
    known = ~np.isnan(field_of(rates, REACTIVITY))

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


def test_the_two_programs_give_the_rate_over_the_control(build, tmp_path):
    treated = build(random_fields(seed=41))
    untreated = build(random_fields(seed=42))
    control = build(random_fields(seed=43))

    output = _normalized(tmp_path, treated, untreated, control)

    treated_rate = field_of(treated, REACTIVITY)
    untreated_rate = field_of(untreated, REACTIVITY)
    control_rate = field_of(control, REACTIVITY)

    with np.errstate(divide="ignore", invalid="ignore"):
        normalized = np.where(control_rate > 0,
                              (treated_rate - untreated_rate) / control_rate,
                              np.float32(np.nan))

    assert np.array_equal(field_of(output, REACTIVITY), normalized, equal_nan=True)


def test_the_two_programs_add_the_errors_in_quadrature(build, tmp_path):
    treated = build(random_fields(seed=44))
    untreated = build(random_fields(seed=45))
    control = build(random_fields(seed=46))

    output = _normalized(tmp_path, treated, untreated, control)

    treated_error = field_of(treated, ERROR)
    untreated_error = field_of(untreated, ERROR)
    control_error = field_of(control, ERROR)
    control_rate = field_of(control, REACTIVITY)
    normalized_rate = field_of(output, REACTIVITY)

    variance = (treated_error ** 2 + untreated_error ** 2
                + (normalized_rate * control_error) ** 2)

    with np.errstate(divide="ignore", invalid="ignore"):
        normalized = np.where(control_rate > 0, np.sqrt(variance) / control_rate,
                              np.float32(np.nan))

    assert np.allclose(field_of(output, ERROR), normalized, rtol=TOLERANCE,
                       equal_nan=True)


def test_every_count_is_summed_over_all_three_inputs(build, tmp_path):
    treated = build({COUNTED: 3}, unmapped=11)
    untreated = build({COUNTED: 5}, unmapped=13)
    control = build({COUNTED: 7}, unmapped=17)

    output = _normalized(tmp_path, treated, untreated, control)

    assert np.all(field_of(output, COUNTED) == 15)
    assert field_of(output, UNMAPPED) == 41


def test_clipping_the_difference_holds_the_normalized_rate_at_zero(build, tmp_path):
    treated = build({REACTIVITY: 0.25, ERROR: 0.1})
    untreated = build({REACTIVITY: 0.75, ERROR: 0.1})
    control = build({REACTIVITY: 0.5, ERROR: 0.1})

    output = _normalized(tmp_path, treated, untreated, control, clip=True)

    assert np.all(field_of(output, REACTIVITY) == 0)


def test_a_clipped_rate_carries_no_uncertainty_from_the_control(build, tmp_path):
    """The error of a ratio scales the control's error by the rate, so at a rate
    of zero the control contributes no term and the error is the quadrature of
    the two inputs divided by the control's rate."""
    treated = build({REACTIVITY: 0.25, ERROR: 0.3})
    untreated = build({REACTIVITY: 0.75, ERROR: 0.4})
    control = build({REACTIVITY: 0.5, ERROR: 0.1})

    output = _normalized(tmp_path, treated, untreated, control, clip=True)

    quadrature = np.sqrt(np.float32(0.3) ** 2 + np.float32(0.4) ** 2)

    assert np.allclose(field_of(output, ERROR), quadrature / np.float32(0.5),
                       rtol=TOLERANCE)
