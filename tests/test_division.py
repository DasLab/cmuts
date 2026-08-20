"""Dividing reactivity rates by a denatured control.

The result depends only on the values in the input files, so the
inputs are written by hand and not counted from an alignment. inputs.py builds
them and outputs.py describes the layout the programs share.
"""

from __future__ import annotations

import numpy as np
import pytest

from combining import DIV_RULES, expected
from inputs import (
    CAP,
    LARGE_COUNT,
    N_REFS,
    NOTES,
    missing_in_each_input,
    not_hdf5,
    random_fields,
    random_values,
)
from outputs import (
    ALL_FIELDS,
    COUNTED,
    COVERAGE,
    ERROR,
    FLOAT_FIELDS,
    REACTIVITY,
    UNMAPPED,
    attributes_of,
    delete_field,
    field_of,
    layout_of,
    read_summary,
    set_field_width,
    write_output,
)
from programs import (
    CMUTS_DIV,
    attempt,
    reported_version,
    run_cmuts,
    run_divide,
    run_subtract,
    try_divide,
)

# The error rounds at every step of the division and the root, so it is compared
# to a relative tolerance. Every other rule rounds once and is exact.
TOLERANCE = 1e-6


@pytest.fixture
def divide(tmp_path):
    """Returns a function that runs the division into a path of its own."""
    def run(rates, control, **options):
        return run_divide(rates, control, tmp_path / "normalized.h5", **options)

    return run


def wanted(name, *inputs):
    """Returns the values cmuts div should write for one field."""
    return expected(DIV_RULES, name, *inputs)


def agrees(result, expectation, name):
    """Returns whether a field matches, allowing a tolerance for the error and
    requiring an exact match for every other field."""
    if name == ERROR:
        return np.allclose(result, expectation, rtol=TOLERANCE, equal_nan=True)

    return np.array_equal(result, expectation)


# ---------------------------------------------------------------------------
# The rules
# ---------------------------------------------------------------------------


def test_the_result_names_the_program_that_wrote_it(build, divide):
    output = divide(build(), build())

    assert attributes_of(output)["program"] == " ".join(CMUTS_DIV)
    assert attributes_of(output)["version"] == reported_version(CMUTS_DIV)


@pytest.mark.parametrize("name", ALL_FIELDS)
def test_each_field_follows_its_rule(build, divide, name):
    rates = build(random_fields(seed=31), unmapped=17)
    control = build(random_fields(seed=32), unmapped=3)

    output = divide(rates, control)

    assert agrees(field_of(output, name), wanted(name, rates, control), name), name


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


def test_a_control_of_zero_leaves_no_reactivity(build, divide):
    output = divide(build({REACTIVITY: 0.5}), build({REACTIVITY: 0.0}))

    assert np.isnan(field_of(output, REACTIVITY)).all()
    assert np.isnan(field_of(output, ERROR)).all()


def test_a_control_below_zero_leaves_no_reactivity(build, divide):
    output = divide(build({REACTIVITY: 0.5}), build({REACTIVITY: -0.25}))

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


def test_counts_stay_whole_and_exact(build, divide):
    rates = build({COUNTED: LARGE_COUNT}, unmapped=LARGE_COUNT)
    control = build({COUNTED: 1}, unmapped=1)

    output = divide(rates, control)

    assert np.all(field_of(output, COUNTED) == LARGE_COUNT + 1)
    assert field_of(output, UNMAPPED) == LARGE_COUNT + 1
    assert field_of(output, COUNTED).dtype == np.dtype("u8")


def test_a_control_is_counted_in_the_totals(build, divide):
    rates = build({COUNTED: 3}, unmapped=11)
    control = build({COUNTED: 5}, unmapped=13)

    output = divide(rates, control)

    assert np.all(field_of(output, COUNTED) == 8)
    assert field_of(output, UNMAPPED) == 24


# ---------------------------------------------------------------------------
# Values that are missing from an input
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", FLOAT_FIELDS)
def test_a_value_either_input_lacks_is_missing_from_the_output(build, divide, name):
    left, right = missing_in_each_input(N_REFS, CAP)

    # The error of a ratio uses both reactivities, so both are set to a known
    # value and the field under test is the only one with a value missing.
    known = {REACTIVITY: 0.5}

    output = divide(build(known | {name: left}), build(known | {name: right}))
    result = field_of(output, name)

    assert np.array_equal(np.isnan(result), np.isnan(left) | np.isnan(right))
    assert not np.isnan(result[0]).any()


def test_the_columns_past_a_reference_stay_nan(build, divide):
    lengths = [6, 4, 2, 1]
    rows = np.full((N_REFS, CAP), np.float32(0.25))

    for row, length in enumerate(lengths):
        rows[row, length:] = np.nan

    output = divide(build({COVERAGE: rows}), build({COVERAGE: rows}))
    result = field_of(output, COVERAGE)

    for row, length in enumerate(lengths):
        assert not np.isnan(result[row, :length]).any(), row
        assert np.isnan(result[row, length:]).all(), row


# ---------------------------------------------------------------------------
# Properties that hold for any two inputs
# ---------------------------------------------------------------------------


def test_the_output_is_shaped_and_typed_like_its_inputs(build, divide):
    rates, control = build(random_fields(seed=2)), build(random_fields(seed=9))

    assert layout_of(divide(rates, control)) == layout_of(rates)


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


def test_two_runs_agree_byte_for_byte(build, tmp_path):
    rates, control = build(random_fields(seed=6)), build(random_fields(seed=60))

    first = run_divide(rates, control, tmp_path / "first.h5")
    second = run_divide(rates, control, tmp_path / "second.h5")

    assert first.read_bytes() == second.read_bytes()


@pytest.mark.parametrize("n_refs, cap", [(1, 1), (1, 40), (3, 1), (400, 2)])
def test_each_field_follows_its_rule_at_any_shape(build, divide, n_refs, cap):
    shaped = dict(n_refs=n_refs, cap=cap)

    rates = build(random_fields(seed=11, **shaped), unmapped=17, **shaped)
    control = build(random_fields(seed=21, **shaped), unmapped=3, **shaped)
    output = divide(rates, control)

    for name in ALL_FIELDS:
        assert agrees(field_of(output, name), wanted(name, rates, control), name), name


# ---------------------------------------------------------------------------
# Inputs that are refused
# ---------------------------------------------------------------------------


def test_inputs_of_different_reference_counts_are_refused(build, tmp_path):
    failed = try_divide(build(n_refs=4), build(n_refs=5), tmp_path / "out.h5")

    assert failed.returncode != 0


def test_inputs_of_different_widths_are_refused(build, tmp_path):
    failed = try_divide(build(cap=6), build(cap=7), tmp_path / "out.h5")

    assert failed.returncode != 0


def test_a_file_holding_no_references_is_refused(tmp_path):
    empty = write_output(tmp_path / "empty.h5", n_refs=0, cap=CAP)

    failed = try_divide(empty, empty, tmp_path / "out.h5")

    assert failed.returncode != 0


def test_something_that_is_not_hdf5_is_refused(build, tmp_path):
    notes = not_hdf5(tmp_path)

    failed = try_divide(notes, build(), tmp_path / "out.h5")

    assert failed.returncode != 0
    assert notes.read_text() == NOTES


# The datasets cmuts div refuses an input for.
REQUIRED = (COVERAGE, REACTIVITY, ERROR)
SKIPPABLE = tuple(name for name in ALL_FIELDS if name not in REQUIRED)


@pytest.mark.parametrize("missing", REQUIRED)
def test_an_input_missing_a_required_dataset_is_refused(build, tmp_path, missing):
    failed = try_divide(delete_field(build(), missing), build(), tmp_path / "out.h5")

    assert failed.returncode != 0

@pytest.mark.parametrize("missing", SKIPPABLE)
def test_an_input_missing_a_dataset_that_is_not_required_is_skipped(build, divide,
                                                                missing):
    output = divide(delete_field(build(), missing), build())

    assert missing not in layout_of(output)
    assert REACTIVITY in layout_of(output)



def test_a_control_missing_a_dataset_is_refused(build, tmp_path):
    failed = try_divide(build(), delete_field(build(), ERROR), tmp_path / "out.h5")

    assert failed.returncode != 0


def test_an_input_whose_datasets_disagree_is_refused(build, tmp_path):
    broken = set_field_width(build(), ERROR, CAP + 3)

    failed = try_divide(broken, build(), tmp_path / "out.h5")

    assert failed.returncode != 0


# ---------------------------------------------------------------------------
# An existing file at the output path
# ---------------------------------------------------------------------------


def test_an_existing_output_is_not_replaced_without_overwrite(build, divide):
    rates, control = build(), build()
    output = divide(rates, control)
    before = output.read_bytes()

    failed = try_divide(rates, control, output)

    assert failed.returncode != 0
    assert output.read_bytes() == before


def test_overwrite_replaces_an_existing_output(build, divide):
    """The two runs are given different inputs, so that the values left at the
    path identify which of them wrote it."""
    output = divide(build({REACTIVITY: 0.25}), build({REACTIVITY: 1.0}))

    rates, control = build({REACTIVITY: 0.75}), build({REACTIVITY: 0.25})
    run_divide(rates, control, output, overwrite=True)

    assert np.all(field_of(output, REACTIVITY) == np.float32(3.0))


@pytest.mark.parametrize("wrong", ["shape", "missing dataset", "not hdf5"])
def test_a_run_that_refuses_its_inputs_leaves_the_output_intact(build, divide,
                                                                tmp_path, wrong):
    bad = {
        "shape": lambda: build(n_refs=N_REFS + 1),
        "missing dataset": lambda: delete_field(build(), REACTIVITY),
        "not hdf5": lambda: not_hdf5(tmp_path),
    }

    output = divide(build(), build())
    before = output.read_bytes()

    failed = try_divide(bad[wrong](), build(), output, overwrite=True)

    assert failed.returncode != 0
    assert output.read_bytes() == before


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------


def test_both_inputs_are_required(build, tmp_path):
    given = attempt([*CMUTS_DIV, "-o", tmp_path / "out.h5", build()])

    assert given.returncode == 2
    assert not (tmp_path / "out.h5").exists()


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


# ---------------------------------------------------------------------------
# End to end
# ---------------------------------------------------------------------------


def test_cmuts_div_reads_what_cmuts_hmm_writes(data, falsifiable, tmp_path):
    """Asserts that the run succeeds and leaves a file shaped like its
    inputs."""
    rates = tmp_path / "rates.h5"
    control = tmp_path / "control.h5"

    summary = read_summary(run_cmuts(data, rates))
    run_cmuts(data, control, min_mapq=30)

    falsifiable(summary.rows > 0)

    output = run_divide(rates, control, tmp_path / "normalized.h5")

    assert layout_of(output) == layout_of(rates)
