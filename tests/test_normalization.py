"""Dividing the rates by one scale taken from the rates themselves.

The scale is pooled over the aggregate rate, one less the product of the three
no-event rates, and divides every rate and error alike. The result depends only on the
values in the input files, so the inputs are written by hand and not counted
from an alignment. inputs.py builds them and outputs.py describes the layout
the programs share. The contracts this program shares with the other readers
of outputs are in test_io.py.

random_fields gives coverage in [0, 1), which no position clears the default
floor with, so every test of the ubr scale sets the coverage it wants.
"""

from __future__ import annotations

import numpy as np
import pytest

from inputs import (
    CAP,
    N_REFS,
    missing_in_each_input,
    not_hdf5,
    random_fields,
)
from normalization import OUTLIER, UBR, expected, factor, pool
from outputs import (
    ALL_FIELDS,
    COUNTED,
    COVERAGE,
    ERROR_FIELDS,
    NORM,
    RATE_FIELDS,
    UNMAPPED,
    field_of,
    layout_of,
)
from programs import CMUTS_NORM, attempt, run_normalize, try_normalize

# The scale is computed in float64 over values narrowed to float32, so a field
# that carries it agrees to a tolerance and not exactly.
TOLERANCE = 1e-6

# Well above the default floor, so every position of an input built with it
# reaches the ubr pool.
COVERED = 1000.0

# The fields the scale divides.
SCALED = RATE_FIELDS + ERROR_FIELDS


@pytest.fixture
def normalize(tmp_path):
    """Returns a function that normalizes any number of inputs together, each
    into an output of its own."""
    def run(*inputs, **options):
        outputs = [tmp_path / f"norm{i}.h5" for i in range(len(inputs))]

        return run_normalize(list(inputs), outputs, **options)

    return run


def covered(values=None, **rest):
    """Field values with the coverage raised above the default floor."""
    return {COVERAGE: COVERED, **(values or {}), **rest}


def every_rate(value):
    """The three channel rates, all set to one value."""
    return dict.fromkeys(RATE_FIELDS, value)


def aggregate(rate: float) -> float:
    """The aggregate of three channels all at one rate."""
    return 1.0 - ((1.0 - rate) ** 3)


def recorded(path) -> float:
    """The scale an output records for itself."""
    return float(field_of(path, NORM))


def besides_the_scale(path) -> dict:
    """The datasets an output holds other than the scale, the one dataset this
    program adds."""
    return {name: shape for name, shape in layout_of(path).items() if name != NORM}


# ---------------------------------------------------------------------------
# The scale
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("scheme", [UBR, OUTLIER])
def test_the_scale_matches_its_oracle(build, normalize, scheme):
    rates = build(covered(random_fields(seed=1)))

    output, = normalize(rates, norm=scheme)

    assert recorded(output) == pytest.approx(factor(scheme, [rates]), rel=TOLERANCE)


def test_a_constant_aggregate_is_its_own_scale(build, normalize):
    """Every value of the pool is the same, so any percentile of it is that
    value: the aggregate of the three rates, not any one of them."""
    output, = normalize(build(covered(every_rate(0.25))))

    assert recorded(output) == pytest.approx(aggregate(0.25), rel=TOLERANCE)

    for name in RATE_FIELDS:
        assert np.allclose(field_of(output, name), 0.25 / aggregate(0.25),
                           rtol=TOLERANCE), name


@pytest.mark.parametrize("scheme", [UBR, OUTLIER])
def test_rates_supporting_no_scale_record_none(build, normalize, scheme):
    """A scale is a divisor, so rates that come to zero support none. The rates are
    left as they are, and the output says it holds no scale rather than one."""
    rates = build(covered(every_rate(0.0)))

    output, = normalize(rates, norm=scheme)

    assert np.isnan(recorded(output))

    for name in RATE_FIELDS:
        assert np.allclose(field_of(output, name), 0.0), name


@pytest.mark.parametrize("scheme", [UBR, OUTLIER])
@pytest.mark.parametrize("name", ALL_FIELDS)
def test_each_field_follows_the_scale(build, normalize, scheme, name):
    rates = build(covered(random_fields(seed=2)), unmapped=17)

    output, = normalize(rates, norm=scheme)
    wanted = expected(rates, name, factor(scheme, [rates]))

    if name in SCALED:
        assert np.allclose(field_of(output, name), wanted, rtol=TOLERANCE,
                           equal_nan=True), name
    else:
        assert np.array_equal(field_of(output, name), wanted), name


def test_every_error_carries_the_scale(build, normalize):
    values = every_rate(0.5) | dict.fromkeys(ERROR_FIELDS, 0.25)

    output, = normalize(build(covered(values)))

    assert recorded(output) == pytest.approx(aggregate(0.5), rel=TOLERANCE)

    for name in ERROR_FIELDS:
        assert np.allclose(field_of(output, name), 0.25 / aggregate(0.5),
                           rtol=TOLERANCE), name


def test_counts_and_coverage_are_left_alone(build, normalize):
    rates = build(covered(random_fields(seed=3)), unmapped=23)

    output, = normalize(rates)

    for name in set(ALL_FIELDS) - set(SCALED):
        assert np.array_equal(field_of(output, name), field_of(rates, name)), name


# ---------------------------------------------------------------------------
# What the pool is drawn from
# ---------------------------------------------------------------------------


def test_a_position_below_the_floor_does_not_reach_the_pool(build, normalize):
    """The floor admits one row and excludes the other, so the scale is the
    admitted row's aggregate alone."""
    coverage = np.full((N_REFS, CAP), np.float32(1.0))
    coverage[0, :] = COVERED

    rates = np.full((N_REFS, CAP), np.float32(0.8), dtype=np.float32)
    rates[0, :] = 0.2

    output, = normalize(build({COVERAGE: coverage} | every_rate(rates)))

    assert recorded(output) == pytest.approx(aggregate(0.2), rel=TOLERANCE)


def test_no_position_clearing_the_floor_leaves_the_rates_alone(build, normalize):
    rates = build(random_fields(seed=4))

    output, = normalize(rates)

    assert recorded(output) == 1.0

    for name in RATE_FIELDS:
        assert np.array_equal(field_of(output, name), field_of(rates, name)), name


def test_lowering_the_floor_admits_more_of_the_pool(build, tmp_path):
    rates = build(random_fields(seed=5))

    strict = run_normalize([rates], [tmp_path / "strict.h5"])
    loose = run_normalize([rates], [tmp_path / "loose.h5"], min_coverage="0")

    assert recorded(strict[0]) == 1.0
    assert recorded(loose[0]) == pytest.approx(factor(UBR, [rates], min_coverage=0),
                                               rel=TOLERANCE)
    assert recorded(loose[0]) != 1.0


def test_the_outlier_scheme_ignores_the_floor(build, tmp_path):
    """Only ubr consults the coverage, so lowering the floor cannot move an
    outlier scale."""
    rates = build(random_fields(seed=6))

    strict = run_normalize([rates], [tmp_path / "strict.h5"], norm=OUTLIER)
    loose = run_normalize([rates], [tmp_path / "loose.h5"], norm=OUTLIER,
                          min_coverage="0")

    assert recorded(strict[0]) == recorded(loose[0])


@pytest.mark.parametrize("channel", RATE_FIELDS)
def test_a_position_missing_any_channel_does_not_reach_the_pool(build, normalize,
                                                                channel):
    """NaN in one channel leaves the position no aggregate, whichever channel
    it is."""
    left, _ = missing_in_each_input(N_REFS, CAP)

    rates = build(covered(every_rate(0.5) | {channel: left}))
    output, = normalize(rates)

    assert pool(UBR, [rates]).size < left.size
    assert recorded(output) == pytest.approx(factor(UBR, [rates]), rel=TOLERANCE)


# ---------------------------------------------------------------------------
# One scale over several inputs
# ---------------------------------------------------------------------------


def test_every_input_is_given_the_same_scale(build, normalize):
    first = build(covered(every_rate(0.2)))
    second = build(covered(every_rate(0.8)))

    outputs = normalize(first, second)

    assert recorded(outputs[0]) == recorded(outputs[1])
    assert recorded(outputs[0]) == pytest.approx(factor(UBR, [first, second]),
                                                 rel=TOLERANCE)


def test_the_pooled_scale_differs_from_either_input_alone(build, normalize, tmp_path):
    """Running the two together is the only way to put them on one scale, which
    is what separate runs give up."""
    first = build(covered(every_rate(0.2)))
    second = build(covered(every_rate(0.8)))

    together = normalize(first, second)
    alone = run_normalize([first], [tmp_path / "alone.h5"])

    assert recorded(alone[0]) == pytest.approx(aggregate(0.2), rel=TOLERANCE)
    assert recorded(together[0]) != pytest.approx(recorded(alone[0]), rel=TOLERANCE)


def test_inputs_of_different_shapes_share_a_scale(build, normalize):
    """The scale is one number, so the inputs need not have been counted against
    the same references."""
    first = build(covered(every_rate(0.4)), n_refs=2, cap=3)
    second = build(covered(every_rate(0.4)), n_refs=5, cap=7)

    outputs = normalize(first, second)

    assert recorded(outputs[0]) == recorded(outputs[1])
    assert besides_the_scale(outputs[0]) == layout_of(first)
    assert besides_the_scale(outputs[1]) == layout_of(second)


def test_each_output_matches_the_input_it_was_paired_with(build, normalize):
    first = build(covered(every_rate(0.2) | {COUNTED: 3}), unmapped=11)
    second = build(covered(every_rate(0.8) | {COUNTED: 5}), unmapped=13)

    outputs = normalize(first, second)

    assert np.all(field_of(outputs[0], COUNTED) == 3)
    assert np.all(field_of(outputs[1], COUNTED) == 5)
    assert field_of(outputs[0], UNMAPPED) == 11
    assert field_of(outputs[1], UNMAPPED) == 13


# ---------------------------------------------------------------------------
# Properties that hold for any input
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n_refs, cap", [(1, 1), (1, 40), (3, 1), (400, 2)])
def test_the_scale_holds_at_any_shape(build, normalize, n_refs, cap):
    shaped = dict(n_refs=n_refs, cap=cap)

    rates = build(covered(random_fields(seed=13, **shaped)), unmapped=7, **shaped)
    output, = normalize(rates)

    assert recorded(output) == pytest.approx(factor(UBR, [rates]), rel=TOLERANCE)

    for name in RATE_FIELDS:
        assert np.allclose(field_of(output, name),
                           expected(rates, name, factor(UBR, [rates])),
                           rtol=TOLERANCE, equal_nan=True), name


# ---------------------------------------------------------------------------
# Inputs that are refused
# ---------------------------------------------------------------------------


def test_fewer_outputs_than_inputs_is_refused(build, tmp_path):
    failed = try_normalize([build(), build()], [tmp_path / "out.h5"])

    assert failed.returncode == 2
    assert not (tmp_path / "out.h5").exists()


def test_more_outputs_than_inputs_is_refused(build, tmp_path):
    failed = try_normalize([build()], [tmp_path / "a.h5", tmp_path / "b.h5"])

    assert failed.returncode == 2
    assert not (tmp_path / "a.h5").exists()


def test_an_input_is_required(tmp_path):
    given = attempt([*CMUTS_NORM, "-o", tmp_path / "out.h5"])

    assert given.returncode == 2
    assert not (tmp_path / "out.h5").exists()


def test_a_later_input_being_refused_writes_no_output(build, tmp_path):
    """Every input is read before the first output is created, so a bad second
    input leaves the first output unwritten."""
    good = build(covered())
    outputs = [tmp_path / "first.h5", tmp_path / "second.h5"]

    failed = try_normalize([good, not_hdf5(tmp_path)], outputs)

    assert failed.returncode != 0
    assert not outputs[0].exists()
    assert not outputs[1].exists()


# ---------------------------------------------------------------------------
# An existing file at the output path
# ---------------------------------------------------------------------------


def test_a_second_output_already_there_leaves_the_first_alone(build, tmp_path):
    """Every output path is checked before any is created."""
    taken = run_normalize([build(covered())], [tmp_path / "taken.h5"])
    fresh = tmp_path / "fresh.h5"

    failed = try_normalize([build(covered()), build(covered())], [fresh, taken[0]])

    assert failed.returncode != 0
    assert not fresh.exists()
