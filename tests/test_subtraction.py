"""Taking an untreated background off a treated run.

The result depends on the values in the input files and on nothing else, so the
inputs are written by hand and not counted from an alignment. outputs.py
describes the layout the two programs share.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from outputs import (
    ALL_FIELDS,
    BY_NAME,
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
from subtraction import expected

# Small enough to write out by hand and to read in a failure, and ragged enough
# that a row, a histogram and a scalar are all of different widths.
N_REFS = 4
CAP = 6

# A count past what a float32 holds exactly, so a value that comes back whole
# was not narrowed to the type the rates use.
LARGE_COUNT = np.uint64(2) ** 40 + 1

NOTES = "months of irreplaceable notes\n"


@pytest.fixture(params=["plain", "chunked"])
def storage(request):
    """The two ways an input may be stored. cmuts-hmm writes chunked, shuffled
    and deflated, and the result must be the same either way, so every test
    that reads values runs against both."""
    return request.param


@pytest.fixture
def build(tmp_path, storage):
    """Returns a function that writes an input file. Each file is named
    separately, so one test may build several."""
    written = itertools.count()

    def make(values=None, *, n_refs=N_REFS, cap=CAP, unmapped=0):
        return write_output(
            tmp_path / f"input{next(written)}.h5",
            n_refs=n_refs, cap=cap, values=values, unmapped=unmapped,
            storage=storage,
        )

    return make


@pytest.fixture
def subtract(tmp_path):
    """Returns a function that runs the subtraction into a path of its own."""
    def run(treated, untreated, **options):
        return run_subtract(treated, untreated, tmp_path / "difference.h5", **options)

    return run


def random_values(field, n_refs=N_REFS, cap=CAP, seed=0):
    """Builds values for one field that differ from column to column and from
    row to row, so that a rule applied along the wrong axis does not agree by
    accident."""
    rng = np.random.default_rng(seed)
    wanted = shape(BY_NAME[field], n_refs, cap)

    if BY_NAME[field].dtype == "u8":
        return rng.integers(0, 1000, size=wanted, dtype=np.uint64)

    return rng.random(size=wanted).astype(np.float32)


def random_fields(seed, n_refs=N_REFS, cap=CAP) -> dict:
    """Builds values for every field at once, so that a test of one field runs
    on a file whose other fields are not zero."""
    return {field.name: random_values(field.name, n_refs, cap, seed=seed + i)
            for i, field in enumerate(FIELDS)}


def missing_in_each_input(rows: int, columns: int):
    """Builds a pair of arrays that are NaN in the treated input only, in the
    untreated input only, and in both."""
    known = np.float32(0.5)

    left = np.full((rows, columns), known, dtype=np.float32)
    right = np.full((rows, columns), known, dtype=np.float32)

    left[1, :] = np.nan
    right[2, :] = np.nan
    left[3, :] = right[3, :] = np.nan

    return left, right


# ---------------------------------------------------------------------------
# The layout, which is shared with cmuts-hmm
# ---------------------------------------------------------------------------


def test_the_layout_written_here_is_the_one_cmuts_hmm_writes(data, falsifiable,
                                                             tmp_path):
    """Checks the description in outputs.py against a real cmuts-hmm run.

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

    assert real[UNMAPPED][0] == (), "a run total belongs to no reference"


def test_the_difference_names_the_program_that_wrote_it(build, subtract):
    output = subtract(build(), build())

    assert attributes_of(output)["program"] == CMUTS_SUB
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
        field_of(output, name), expected(name, treated, untreated),
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

    assert (unclipped < 0).any(), "nothing was negative, so nothing was tested"
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
# The denatured control
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ALL_FIELDS)
def test_each_field_follows_its_rule_against_a_control(build, subtract, name):
    """The error is compared to a tolerance because it rounds eight times over
    a division and a root; every other rule rounds once and is exact."""
    treated = build(random_fields(seed=31), unmapped=17)
    untreated = build(random_fields(seed=32), unmapped=3)
    denatured = build(random_fields(seed=33), unmapped=5)

    output = subtract(treated, untreated, denatured=denatured)

    result = field_of(output, name)
    wanted = expected(name, treated, untreated, denatured)

    if name == ERROR:
        assert np.allclose(result, wanted, rtol=1e-6, equal_nan=True), name
    else:
        assert np.array_equal(result, wanted), name


def test_a_control_divides_the_difference(build, subtract):
    treated = build({REACTIVITY: 0.75})
    untreated = build({REACTIVITY: 0.25})
    denatured = build({REACTIVITY: 0.5})

    output = subtract(treated, untreated, denatured=denatured)

    assert np.all(field_of(output, REACTIVITY) == np.float32(1.0))


def test_a_control_of_zero_leaves_no_reactivity(build, subtract):
    treated = build({REACTIVITY: 0.75})
    untreated = build({REACTIVITY: 0.25})
    denatured = build({REACTIVITY: 0.0})

    output = subtract(treated, untreated, denatured=denatured)

    assert np.isnan(field_of(output, REACTIVITY)).all()
    assert np.isnan(field_of(output, ERROR)).all()


def test_a_control_of_nan_leaves_no_reactivity(build, subtract):
    control = np.full((N_REFS, CAP), np.float32(0.5), dtype=np.float32)
    control[2, :] = np.nan

    output = subtract(build({REACTIVITY: 0.75}), build({REACTIVITY: 0.25}),
                      denatured=build({REACTIVITY: control}))
    result = field_of(output, REACTIVITY)

    assert np.isnan(result[2]).all()
    assert not np.isnan(result[0]).any()


def test_a_control_of_ones_leaves_the_difference_alone(build, tmp_path):
    treated, untreated = build(random_fields(seed=34)), build(random_fields(seed=35))
    ones = build({REACTIVITY: 1.0})

    plain = run_subtract(treated, untreated, tmp_path / "plain.h5")
    controlled = run_subtract(treated, untreated, tmp_path / "controlled.h5",
                              denatured=ones)

    assert np.array_equal(field_of(controlled, REACTIVITY),
                          field_of(plain, REACTIVITY))


def test_a_control_can_only_widen_the_error(build, tmp_path):
    treated, untreated = build(random_fields(seed=36)), build(random_fields(seed=37))
    ones = build({REACTIVITY: 1.0, ERROR: random_values(ERROR, seed=38)})

    plain = run_subtract(treated, untreated, tmp_path / "plain.h5")
    controlled = run_subtract(treated, untreated, tmp_path / "controlled.h5",
                              denatured=ones)

    assert np.all(field_of(controlled, ERROR) >= field_of(plain, ERROR))


def test_a_control_is_counted_in_the_totals(build, subtract):
    treated = build({COUNTED: 3}, unmapped=11)
    untreated = build({COUNTED: 5}, unmapped=13)
    denatured = build({COUNTED: 7}, unmapped=17)

    output = subtract(treated, untreated, denatured=denatured)

    assert np.all(field_of(output, COUNTED) == 15)
    assert field_of(output, UNMAPPED) == 41


def test_clipping_holds_a_ratio_at_zero(build, subtract):
    treated = build({REACTIVITY: 0.25})
    untreated = build({REACTIVITY: 0.75})
    denatured = build({REACTIVITY: 0.5})

    output = subtract(treated, untreated, denatured=denatured, clip=True)

    assert np.all(field_of(output, REACTIVITY) == 0)


# The ways a control can disagree with the inputs, and the word each is
# reported with.
DISAGREEING_CONTROLS = {
    "references": lambda build: build(n_refs=N_REFS + 1),
    "wide": lambda build: build(cap=CAP + 1),
}


@pytest.mark.parametrize("wrong", sorted(DISAGREEING_CONTROLS))
def test_a_control_disagreeing_with_the_inputs_is_refused(build, tmp_path, wrong):
    control = DISAGREEING_CONTROLS[wrong](build)

    failed = try_subtract(build(), build(), tmp_path / "out.h5", denatured=control)

    assert failed.returncode != 0


def test_a_control_missing_a_dataset_is_refused(build, tmp_path):
    failed = try_subtract(build(), build(), tmp_path / "out.h5",
                          denatured=delete_field(build(), ERROR))

    assert failed.returncode != 0


# ---------------------------------------------------------------------------
# Values that were never measured
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", FLOAT_FIELDS)
def test_a_value_either_input_lacks_is_missing_from_the_output(build, subtract, name):
    left, right = missing_in_each_input(N_REFS, CAP)

    output = subtract(build({name: left}), build({name: right}))
    result = field_of(output, name)

    assert np.array_equal(np.isnan(result), np.isnan(left) | np.isnan(right))
    assert not np.isnan(result[0]).any(), "nothing is lost where both are known"


def test_the_columns_past_a_reference_stay_nan(build, subtract):
    lengths = [6, 4, 2, 1]
    rows = np.full((N_REFS, CAP), np.float32(0.25))

    for row, length in enumerate(lengths):
        rows[row, length:] = np.nan

    output = subtract(build({COVERAGE: rows}), build({COVERAGE: rows}))
    result = field_of(output, COVERAGE)

    for row, length in enumerate(lengths):
        assert not np.isnan(result[row, :length]).any(), f"row {row} within"
        assert np.isnan(result[row, length:]).all(), f"row {row} past its end"


# ---------------------------------------------------------------------------
# What holds of any two files
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

    assert known.any(), "nothing was known, so nothing was tested"
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
    nothing = build({REACTIVITY: 0.0, ERROR: 0.0})

    output = subtract(treated, nothing)

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
            field_of(output, name), expected(name, treated, untreated),
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
    notes = tmp_path / "notes.txt"
    notes.write_text(NOTES)

    failed = try_subtract(notes, build(), tmp_path / "out.h5")

    assert failed.returncode != 0
    assert notes.read_text() == NOTES


@pytest.mark.parametrize("missing", ALL_FIELDS)
def test_an_input_missing_any_dataset_is_refused(build, tmp_path, missing):
    failed = try_subtract(delete_field(build(), missing), build(), tmp_path / "out.h5")

    assert failed.returncode != 0


def test_an_input_whose_datasets_disagree_is_refused(build, tmp_path):
    broken = set_field_width(build(), ERROR, CAP + 3)

    failed = try_subtract(broken, build(), tmp_path / "out.h5")

    assert failed.returncode != 0


# ---------------------------------------------------------------------------
# What is already at the output path
# ---------------------------------------------------------------------------


def test_an_existing_output_is_not_replaced_without_overwrite(build, subtract):
    treated, untreated = build(), build()
    output = subtract(treated, untreated)
    before = output.read_bytes()

    failed = try_subtract(treated, untreated, output)

    assert failed.returncode != 0
    assert output.read_bytes() == before, "the first result is untouched"


def test_overwrite_replaces_an_existing_output(build, subtract):
    """The two runs are given different inputs, so that the values left at the
    path identify which of them wrote it."""
    output = subtract(build({REACTIVITY: 0.25}), build({REACTIVITY: 0.0}))

    treated, untreated = build({REACTIVITY: 0.75}), build({REACTIVITY: 0.25})
    run_subtract(treated, untreated, output, overwrite=True)

    assert np.all(field_of(output, REACTIVITY) == np.float32(0.5))


def _not_hdf5(tmp_path):
    notes = tmp_path / "notes.txt"
    notes.write_text(NOTES)

    return notes


# The ways an input can be refused, each built from the build fixture and the
# temporary directory.
BAD_INPUTS = {
    "shape": lambda build, tmp_path: build(n_refs=N_REFS + 1),
    "missing dataset": lambda build, tmp_path: delete_field(build(), UNMAPPED),
    "not hdf5": lambda build, tmp_path: _not_hdf5(tmp_path),
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
    given = attempt([CMUTS_SUB, "-o", tmp_path / "out.h5", build()])

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

    summary = run_cmuts(data, treated)
    run_cmuts(data, untreated, min_mapq=30)

    falsifiable(summary.rows > 0)

    output = run_subtract(treated, untreated, tmp_path / "difference.h5")

    assert layout_of(output) == layout_of(treated)
