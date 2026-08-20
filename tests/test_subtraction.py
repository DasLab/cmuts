"""Taking an untreated background off a treated run.

The result depends only on the values in the input files, so the
inputs are written by hand and not counted from an alignment. inputs.py builds
them and outputs.py describes the layout the programs share.
"""

from __future__ import annotations

import numpy as np
import pytest

from combining import SUB_RULES, expected
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
    FIELDS,
    FLOAT_FIELDS,
    REACTIVITY,
    UNMAPPED,
    attributes_of,
    delete_field,
    field_of,
    layout_of,
    read_summary,
    set_field_width,
    shape,
    write_output,
)
from programs import (
    CMUTS_SUB,
    attempt,
    reported_version,
    run_cmuts,
    run_subtract,
    try_subtract,
)


@pytest.fixture
def subtract(tmp_path):
    """Returns a function that runs the subtraction into a path of its own."""
    def run(treated, untreated, **options):
        return run_subtract(treated, untreated, tmp_path / "difference.h5", **options)

    return run


def wanted(name, *inputs):
    """Returns the values cmuts sub should write for one field."""
    return expected(SUB_RULES, name, *inputs)


# ---------------------------------------------------------------------------
# The layout, which is shared with cmuts hmm
# ---------------------------------------------------------------------------


def test_the_layout_written_here_is_the_one_cmuts_hmm_writes(data, falsifiable,
                                                             tmp_path):
    """Checks the description in outputs.py against a real cmuts hmm run.

    Compares names, types and widths, which keeps the reactivity calculation
    out of the comparison.
    """
    counted = tmp_path / "counted.h5"
    run_cmuts(data, counted)
    real = layout_of(counted)

    assert set(real) == set(ALL_FIELDS)

    n_refs, cap = real[COVERAGE][0]

    # The widths compared below come from the rows, one per reference.
    falsifiable(n_refs > 0)

    for field in FIELDS:
        found, dtype = real[field.name]

        assert found == shape(field, n_refs, cap), field.name
        assert dtype == np.dtype(field.dtype), field.name

    assert real[UNMAPPED][0] == ()


def test_the_difference_names_the_program_that_wrote_it(build, subtract):
    output = subtract(build(), build())

    assert attributes_of(output)["program"] == " ".join(CMUTS_SUB)
    assert attributes_of(output)["version"] == reported_version(CMUTS_SUB)


# ---------------------------------------------------------------------------
# The rules
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ALL_FIELDS)
def test_each_field_follows_its_rule(build, subtract, name):
    treated = build(random_fields(seed=1), unmapped=17)
    untreated = build(random_fields(seed=50), unmapped=3)

    output = subtract(treated, untreated)

    assert np.array_equal(
        field_of(output, name), wanted(name, treated, untreated),
    ), name


def test_a_background_above_the_signal_leaves_a_negative_reactivity(build, subtract):
    treated = build({REACTIVITY: 0.25})
    untreated = build({REACTIVITY: 0.75})

    difference = field_of(subtract(treated, untreated), REACTIVITY)

    assert np.all(difference == np.float32(-0.5))


def test_clipping_holds_the_difference_at_zero(build, subtract):
    treated = build({REACTIVITY: 0.25})
    untreated = build({REACTIVITY: 0.75})

    output = subtract(treated, untreated, clip=True)

    assert np.all(field_of(output, REACTIVITY) == 0)


def test_clipping_leaves_a_difference_above_zero_alone(build, tmp_path):
    treated, untreated = build(random_fields(seed=13)), build(random_fields(seed=14))

    plain = run_subtract(treated, untreated, tmp_path / "plain.h5")
    clipped = run_subtract(treated, untreated, tmp_path / "clipped.h5", clip=True)

    unclipped = field_of(plain, REACTIVITY)

    assert (unclipped < 0).any()
    assert np.array_equal(field_of(clipped, REACTIVITY), np.maximum(unclipped, 0))


def test_clipping_does_not_raise_a_missing_value_to_zero(build, subtract):
    left, right = missing_in_each_input(N_REFS, CAP)

    output = subtract(build({REACTIVITY: left}), build({REACTIVITY: right}), clip=True)
    result = field_of(output, REACTIVITY)

    assert np.array_equal(np.isnan(result), np.isnan(left) | np.isnan(right))


def test_clipping_reaches_no_field_but_the_reactivity(build, tmp_path):
    treated = build(random_fields(seed=15), unmapped=11)
    untreated = build(random_fields(seed=16), unmapped=4)

    plain = run_subtract(treated, untreated, tmp_path / "plain.h5")
    clipped = run_subtract(treated, untreated, tmp_path / "clipped.h5", clip=True)

    for name in set(ALL_FIELDS) - {REACTIVITY}:
        assert np.array_equal(field_of(clipped, name), field_of(plain, name)), name


def test_an_error_is_never_reduced_by_subtracting(build, subtract):
    treated = build({ERROR: random_values(ERROR, seed=7)})
    untreated = build({ERROR: random_values(ERROR, seed=8)})

    output = subtract(treated, untreated)

    assert np.all(field_of(output, ERROR) >= field_of(treated, ERROR))
    assert np.all(field_of(output, ERROR) >= field_of(untreated, ERROR))


def test_counts_stay_whole_and_exact(build, subtract):
    treated = build({COUNTED: LARGE_COUNT}, unmapped=LARGE_COUNT)
    untreated = build({COUNTED: 1}, unmapped=1)

    output = subtract(treated, untreated)

    assert np.all(field_of(output, COUNTED) == LARGE_COUNT + 1)
    assert field_of(output, UNMAPPED) == LARGE_COUNT + 1
    assert field_of(output, COUNTED).dtype == np.dtype("u8")


# ---------------------------------------------------------------------------
# Values that are missing from an input
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", FLOAT_FIELDS)
def test_a_value_either_input_lacks_is_missing_from_the_output(build, subtract, name):
    left, right = missing_in_each_input(N_REFS, CAP)

    output = subtract(build({name: left}), build({name: right}))
    result = field_of(output, name)

    assert np.array_equal(np.isnan(result), np.isnan(left) | np.isnan(right))
    assert not np.isnan(result[0]).any()


def test_the_columns_past_a_reference_stay_nan(build, subtract):
    lengths = [6, 4, 2, 1]
    rows = np.full((N_REFS, CAP), np.float32(0.25))

    for row, length in enumerate(lengths):
        rows[row, length:] = np.nan

    output = subtract(build({COVERAGE: rows}), build({COVERAGE: rows}))
    result = field_of(output, COVERAGE)

    for row, length in enumerate(lengths):
        assert not np.isnan(result[row, :length]).any(), row
        assert np.isnan(result[row, length:]).all(), row


# ---------------------------------------------------------------------------
# Properties that hold for any two inputs
# ---------------------------------------------------------------------------


def test_the_output_is_shaped_and_typed_like_its_inputs(build, subtract):
    treated, untreated = build(random_fields(seed=2)), build(random_fields(seed=9))

    assert layout_of(subtract(treated, untreated)) == layout_of(treated)


def test_a_file_against_itself_leaves_a_reactivity_of_zero(build, subtract):
    values = random_fields(seed=3)
    values[REACTIVITY][2, 3] = np.nan

    treated = build(values)
    result = field_of(subtract(treated, treated), REACTIVITY)
    known = ~np.isnan(field_of(treated, REACTIVITY))

    assert known.any()
    assert np.all(result[known] == 0)
    assert np.isnan(result[~known]).all()


def test_swapping_the_two_inputs_negates_only_the_reactivity(build, subtract, tmp_path):
    treated, untreated = build(random_fields(seed=4)), build(random_fields(seed=40))

    forward = subtract(treated, untreated)
    backward = run_subtract(untreated, treated, tmp_path / "backward.h5")

    assert np.array_equal(field_of(backward, REACTIVITY),
                          -field_of(forward, REACTIVITY))

    for name in set(ALL_FIELDS) - {REACTIVITY}:
        assert np.array_equal(field_of(backward, name), field_of(forward, name)), name


def test_a_background_of_zeros_leaves_the_treated_run_unchanged(build, subtract):
    treated = build(random_fields(seed=5))
    zeros = build({REACTIVITY: 0.0, ERROR: 0.0})

    output = subtract(treated, zeros)

    for name in ALL_FIELDS:
        assert np.array_equal(field_of(output, name), field_of(treated, name)), name


def test_two_runs_agree_byte_for_byte(build, tmp_path):
    treated, untreated = build(random_fields(seed=6)), build(random_fields(seed=60))

    first = run_subtract(treated, untreated, tmp_path / "first.h5")
    second = run_subtract(treated, untreated, tmp_path / "second.h5")

    assert first.read_bytes() == second.read_bytes()


@pytest.mark.parametrize("n_refs, cap", [(1, 1), (1, 40), (3, 1), (400, 2)])
def test_each_field_follows_its_rule_at_any_shape(build, subtract, n_refs, cap):
    shaped = dict(n_refs=n_refs, cap=cap)

    treated = build(random_fields(seed=11, **shaped), unmapped=17, **shaped)
    untreated = build(random_fields(seed=21, **shaped), unmapped=3, **shaped)
    output = subtract(treated, untreated)

    for name in ALL_FIELDS:
        assert np.array_equal(
            field_of(output, name), wanted(name, treated, untreated),
        ), name


# ---------------------------------------------------------------------------
# Inputs that are refused
# ---------------------------------------------------------------------------


def test_inputs_of_different_reference_counts_are_refused(build, tmp_path):
    failed = try_subtract(build(n_refs=4), build(n_refs=5), tmp_path / "out.h5")

    assert failed.returncode != 0


def test_inputs_of_different_widths_are_refused(build, tmp_path):
    failed = try_subtract(build(cap=6), build(cap=7), tmp_path / "out.h5")

    assert failed.returncode != 0


def test_a_file_holding_no_references_is_refused(tmp_path):
    empty = write_output(tmp_path / "empty.h5", n_refs=0, cap=CAP)

    failed = try_subtract(empty, empty, tmp_path / "out.h5")

    assert failed.returncode != 0


def test_something_that_is_not_hdf5_is_refused(build, tmp_path):
    notes = not_hdf5(tmp_path)

    failed = try_subtract(notes, build(), tmp_path / "out.h5")

    assert failed.returncode != 0
    assert notes.read_text() == NOTES


# The datasets cmuts sub refuses an input for.
REQUIRED = (COVERAGE, REACTIVITY)
SKIPPABLE = tuple(name for name in ALL_FIELDS if name not in REQUIRED)


@pytest.mark.parametrize("missing", REQUIRED)
def test_an_input_missing_a_required_dataset_is_refused(build, tmp_path, missing):
    failed = try_subtract(delete_field(build(), missing), build(), tmp_path / "out.h5")

    assert failed.returncode != 0

@pytest.mark.parametrize("missing", SKIPPABLE)
def test_an_input_missing_a_dataset_that_is_not_required_is_skipped(build, subtract,
                                                                missing):
    output = subtract(delete_field(build(), missing), build())

    assert missing not in layout_of(output)
    assert REACTIVITY in layout_of(output)



def test_an_input_whose_datasets_disagree_is_refused(build, tmp_path):
    broken = set_field_width(build(), ERROR, CAP + 3)

    failed = try_subtract(broken, build(), tmp_path / "out.h5")

    assert failed.returncode != 0


# ---------------------------------------------------------------------------
# An existing file at the output path
# ---------------------------------------------------------------------------


def test_an_existing_output_is_not_replaced_without_overwrite(build, subtract):
    treated, untreated = build(), build()
    output = subtract(treated, untreated)
    before = output.read_bytes()

    failed = try_subtract(treated, untreated, output)

    assert failed.returncode != 0
    assert output.read_bytes() == before


def test_overwrite_replaces_an_existing_output(build, subtract):
    """The two runs are given different inputs, so that the values left at the
    path identify which of them wrote it."""
    output = subtract(build({REACTIVITY: 0.25}), build({REACTIVITY: 0.0}))

    treated, untreated = build({REACTIVITY: 0.75}), build({REACTIVITY: 0.25})
    run_subtract(treated, untreated, output, overwrite=True)

    assert np.all(field_of(output, REACTIVITY) == np.float32(0.5))


# The ways an input can be refused, each built from the build fixture and the
# temporary directory.
BAD_INPUTS = {
    "shape": lambda build, tmp_path: build(n_refs=N_REFS + 1),
    "missing dataset": lambda build, tmp_path: delete_field(build(), REACTIVITY),
    "not hdf5": lambda build, tmp_path: not_hdf5(tmp_path),
}


@pytest.mark.parametrize("wrong", sorted(BAD_INPUTS))
def test_a_run_that_refuses_its_inputs_leaves_the_output_intact(build, subtract,
                                                                tmp_path, wrong):
    output = subtract(build(), build())
    before = output.read_bytes()

    bad = BAD_INPUTS[wrong](build, tmp_path)

    failed = try_subtract(bad, build(), output, overwrite=True)

    assert failed.returncode != 0
    assert output.read_bytes() == before


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------


def test_both_inputs_are_required(build, tmp_path):
    given = attempt([*CMUTS_SUB, "-o", tmp_path / "out.h5", build()])

    assert given.returncode == 2
    assert not (tmp_path / "out.h5").exists()


# ---------------------------------------------------------------------------
# End to end
# ---------------------------------------------------------------------------


def test_cmuts_sub_reads_what_cmuts_hmm_writes(data, falsifiable, tmp_path):
    """Asserts that the run succeeds and leaves a file shaped like its
    inputs."""
    treated = tmp_path / "treated.h5"
    untreated = tmp_path / "untreated.h5"

    summary = read_summary(run_cmuts(data, treated))
    run_cmuts(data, untreated, min_mapq=30)

    falsifiable(summary.rows > 0)

    output = run_subtract(treated, untreated, tmp_path / "difference.h5")

    assert layout_of(output) == layout_of(treated)
