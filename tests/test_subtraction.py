"""Taking an untreated background off a treated run.

The result depends on the values in the input files and on nothing else, so the
inputs are written by hand rather than counted from an alignment. outputs.py
describes the layout the two programs share.
"""

from __future__ import annotations

import itertools
import subprocess

import numpy as np
import pytest

from outputs import (
    BY_NAME,
    FIELDS,
    RUN_TOTAL,
    combined,
    datasets_of,
    field_of,
    rewidened,
    shape,
    without,
    write_output,
)
from support import CMUTS_SUB, run_cmuts, run_subtract, try_subtract

# Small enough to write out by hand and to read in a failure, and ragged enough
# that a row, a histogram and a scalar are all of different widths.
N_REFS = 4
CAP = 6

# Every dataset the arithmetic applies to, run totals included.
COMBINED = [field.name for field in FIELDS] + [RUN_TOTAL]

# The datasets holding a rate or a fraction, which are the ones that can be
# marked as never measured.
FLOATING = ("coverage", "reactivity", "error")


@pytest.fixture(params=["plain", "chunked"])
def storage(request):
    """The two ways an input may be stored. cmuts writes chunked, shuffled and
    deflated, and the result must be the same either way, so every test that
    reads values runs against both."""
    return request.param


@pytest.fixture
def build(tmp_path, storage):
    """Writes an input, returning its path. Each is named separately, so one
    test may build several."""
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
    """Runs the subtraction into a path of its own."""
    def run(treated, untreated, **options):
        return run_subtract(treated, untreated, tmp_path / "difference.h5", **options)

    return run


def spread(field, n_refs=N_REFS, cap=CAP, seed=0):
    """Values filling one field, differing from column to column and from row
    to row, so that a rule applied along the wrong axis does not agree by
    accident."""
    rng = np.random.default_rng(seed)
    wanted = shape(BY_NAME[field], n_refs, cap)

    if BY_NAME[field].dtype == "u8":
        return rng.integers(0, 1000, size=wanted, dtype=np.uint64)

    return rng.random(size=wanted).astype(np.float32)


def everything(seed):
    """A value for every field at once, so that a test of one field runs on a
    file whose other fields are not zero."""
    return {field.name: spread(field.name, seed=seed + i)
            for i, field in enumerate(FIELDS)}


# ---------------------------------------------------------------------------
# The layout, which is shared with cmuts
# ---------------------------------------------------------------------------


def test_the_layout_written_here_is_the_one_cmuts_writes(data, tmp_path):
    """Checks the description in outputs.py against a real cmuts run.

    Compares names, types and widths and never a value, which keeps the
    reactivity calculation out of the comparison.
    """
    counted = tmp_path / "counted.h5"
    run_cmuts(data, counted)
    real = datasets_of(counted)

    assert set(real) == {field.name for field in FIELDS} | {RUN_TOTAL}

    n_refs, cap = real["coverage"][0]

    for field in FIELDS:
        found, dtype = real[field.name]

        assert found == shape(field, n_refs, cap), field.name
        assert dtype == np.dtype(field.dtype), field.name

    assert real[RUN_TOTAL][0] == (), "a run total belongs to no reference"


# ---------------------------------------------------------------------------
# The rules
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", COMBINED)
def test_each_field_follows_its_rule(build, subtract, name):
    treated = build(everything(seed=1), unmapped=17)
    untreated = build(everything(seed=50), unmapped=3)

    output = subtract(treated, untreated)

    assert np.array_equal(
        field_of(output, name), combined(name, treated, untreated),
    ), name


def test_a_background_above_the_signal_leaves_a_negative_reactivity(build, subtract):
    treated = build({"reactivity": 0.25})
    untreated = build({"reactivity": 0.75})

    difference = field_of(subtract(treated, untreated), "reactivity")

    assert np.all(difference == np.float32(-0.5))


def test_clipping_holds_the_difference_at_zero(build, subtract):
    treated = build({"reactivity": 0.25})
    untreated = build({"reactivity": 0.75})

    assert np.all(field_of(subtract(treated, untreated, clip=True), "reactivity") == 0)


def test_clipping_leaves_a_difference_above_zero_alone(build, tmp_path):
    """Compared against the same pair of inputs run without the flag."""
    treated, untreated = build(everything(seed=13)), build(everything(seed=14))

    plain = run_subtract(treated, untreated, tmp_path / "plain.h5")
    clipped = run_subtract(treated, untreated, tmp_path / "clipped.h5", clip=True)

    unclipped = field_of(plain, "reactivity")

    assert (unclipped < 0).any(), "nothing was negative, so nothing was tested"
    assert np.array_equal(field_of(clipped, "reactivity"), np.maximum(unclipped, 0))


def test_clipping_does_not_raise_a_missing_value_to_zero(build, subtract):
    known, missing = np.float32(0.25), np.float32(np.nan)

    left = np.full((N_REFS, CAP), known, dtype=np.float32)
    right = np.full((N_REFS, CAP), known, dtype=np.float32)

    left[1, :] = missing
    right[2, :] = missing
    left[3, :] = right[3, :] = missing

    output = subtract(build({"reactivity": left}), build({"reactivity": right}),
                      clip=True)
    result = field_of(output, "reactivity")

    assert np.array_equal(np.isnan(result), np.isnan(left) | np.isnan(right))


def test_clipping_reaches_no_field_but_the_reactivity(build, tmp_path):
    treated = build(everything(seed=15), unmapped=11)
    untreated = build(everything(seed=16), unmapped=4)

    plain = run_subtract(treated, untreated, tmp_path / "plain.h5")
    clipped = run_subtract(treated, untreated, tmp_path / "clipped.h5", clip=True)

    for name in set(COMBINED) - {"reactivity"}:
        assert np.array_equal(field_of(clipped, name), field_of(plain, name)), name


def test_an_error_is_never_reduced_by_subtracting(build, subtract):
    treated = build({"error": spread("error", seed=7)})
    untreated = build({"error": spread("error", seed=8)})

    output = subtract(treated, untreated)

    assert np.all(field_of(output, "error") >= field_of(treated, "error"))
    assert np.all(field_of(output, "error") >= field_of(untreated, "error"))


def test_counts_stay_whole_and_exact(build, subtract):
    """A count of 2**40 is past what a float32 holds exactly, so coming back
    whole confirms it was not narrowed to the type the rates use."""
    large = np.uint64(2) ** 40 + 1

    treated = build({"reads/counted": large}, unmapped=large)
    untreated = build({"reads/counted": 1}, unmapped=1)

    output = subtract(treated, untreated)

    assert np.all(field_of(output, "reads/counted") == large + 1)
    assert field_of(output, RUN_TOTAL) == large + 1
    assert field_of(output, "reads/counted").dtype == np.dtype("u8")


# ---------------------------------------------------------------------------
# The denatured control
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", COMBINED)
def test_each_field_follows_its_rule_against_a_control(build, subtract, name):
    """The error is compared to a tolerance because it rounds eight times over
    a division and a root; every other rule rounds once and is exact."""
    treated = build(everything(seed=31), unmapped=17)
    untreated = build(everything(seed=32), unmapped=3)
    denatured = build(everything(seed=33), unmapped=5)

    output = subtract(treated, untreated, denatured=denatured)

    result = field_of(output, name)
    wanted = combined(name, treated, untreated, denatured)

    if name == "error":
        assert np.allclose(result, wanted, rtol=1e-6, equal_nan=True), name
    else:
        assert np.array_equal(result, wanted), name


def test_a_control_divides_the_difference(build, subtract):
    treated = build({"reactivity": 0.75})
    untreated = build({"reactivity": 0.25})
    denatured = build({"reactivity": 0.5})

    output = subtract(treated, untreated, denatured=denatured)

    assert np.all(field_of(output, "reactivity") == np.float32(1.0))


def test_a_control_that_measured_nothing_leaves_no_reactivity(build, subtract):
    treated = build({"reactivity": 0.75})
    untreated = build({"reactivity": 0.25})
    denatured = build({"reactivity": 0.0})

    output = subtract(treated, untreated, denatured=denatured)

    assert np.isnan(field_of(output, "reactivity")).all()
    assert np.isnan(field_of(output, "error")).all()


def test_a_control_of_nan_leaves_no_reactivity(build, subtract):
    known, missing = np.float32(0.5), np.float32(np.nan)

    control = np.full((N_REFS, CAP), known, dtype=np.float32)
    control[2, :] = missing

    output = subtract(build({"reactivity": 0.75}), build({"reactivity": 0.25}),
                      denatured=build({"reactivity": control}))
    result = field_of(output, "reactivity")

    assert np.isnan(result[2]).all()
    assert not np.isnan(result[0]).any()


def test_a_control_of_ones_leaves_the_difference_alone(build, tmp_path):
    treated, untreated = build(everything(seed=34)), build(everything(seed=35))
    ones = build({"reactivity": 1.0})

    plain = run_subtract(treated, untreated, tmp_path / "plain.h5")
    controlled = run_subtract(treated, untreated, tmp_path / "controlled.h5",
                              denatured=ones)

    assert np.array_equal(field_of(controlled, "reactivity"),
                          field_of(plain, "reactivity"))


def test_a_control_can_only_widen_the_error(build, subtract, tmp_path):
    treated, untreated = build(everything(seed=36)), build(everything(seed=37))
    ones = build({"reactivity": 1.0, "error": spread("error", seed=38)})

    plain = run_subtract(treated, untreated, tmp_path / "plain.h5")
    controlled = run_subtract(treated, untreated, tmp_path / "controlled.h5",
                              denatured=ones)

    assert np.all(field_of(controlled, "error") >= field_of(plain, "error"))


def test_a_control_is_counted_in_the_totals(build, subtract):
    treated = build({"reads/counted": 3}, unmapped=11)
    untreated = build({"reads/counted": 5}, unmapped=13)
    denatured = build({"reads/counted": 7}, unmapped=17)

    output = subtract(treated, untreated, denatured=denatured)

    assert np.all(field_of(output, "reads/counted") == 15)
    assert field_of(output, RUN_TOTAL) == 41


def test_clipping_holds_a_ratio_at_zero(build, subtract):
    treated = build({"reactivity": 0.25})
    untreated = build({"reactivity": 0.75})
    denatured = build({"reactivity": 0.5})

    output = subtract(treated, untreated, denatured=denatured, clip=True)

    assert np.all(field_of(output, "reactivity") == 0)


@pytest.mark.parametrize("wrong", ["references", "wide"])
def test_a_control_disagreeing_with_the_inputs_is_refused(build, tmp_path, wrong):
    control = {"references": lambda: build(n_refs=N_REFS + 1),
               "wide": lambda: build(cap=CAP + 1)}[wrong]()

    attempt = try_subtract(build(), build(), tmp_path / "out.h5",
                           denatured=control)

    assert attempt.returncode != 0
    assert wrong in attempt.stderr


def test_a_control_missing_a_dataset_is_refused(build, tmp_path):
    attempt = try_subtract(build(), build(), tmp_path / "out.h5",
                           denatured=without(build(), "error"))

    assert attempt.returncode != 0
    assert "error" in attempt.stderr


# ---------------------------------------------------------------------------
# What was never measured
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", FLOATING)
def test_a_value_either_input_lacks_is_missing_from_the_output(build, subtract, name):
    known, missing = np.float32(0.5), np.float32(np.nan)

    left = np.full((N_REFS, CAP), known, dtype=np.float32)
    right = np.full((N_REFS, CAP), known, dtype=np.float32)

    left[1, :] = missing                    # missing in the treated run alone
    right[2, :] = missing                   # missing in the background alone
    left[3, :] = right[3, :] = missing      # missing in both

    output = subtract(build({name: left}), build({name: right}))
    result = field_of(output, name)

    assert np.array_equal(np.isnan(result), np.isnan(left) | np.isnan(right))
    assert not np.isnan(result[0]).any(), "nothing is lost where both are known"


def test_the_columns_past_a_reference_stay_nan(build, subtract):
    lengths = [6, 4, 2, 1]
    rows = np.full((N_REFS, CAP), np.float32(0.25))

    for row, length in enumerate(lengths):
        rows[row, length:] = np.nan

    output = subtract(build({"coverage": rows}), build({"coverage": rows}))
    result = field_of(output, "coverage")

    for row, length in enumerate(lengths):
        assert not np.isnan(result[row, :length]).any(), f"row {row} within"
        assert np.isnan(result[row, length:]).all(), f"row {row} past its end"


# ---------------------------------------------------------------------------
# What holds of any two files
# ---------------------------------------------------------------------------


def test_the_output_is_shaped_and_typed_like_its_inputs(build, subtract):
    treated, untreated = build(everything(seed=2)), build(everything(seed=9))

    assert datasets_of(subtract(treated, untreated)) == datasets_of(treated)


def test_a_file_against_itself_leaves_a_reactivity_of_zero(build, subtract):
    values = everything(seed=3)
    values["reactivity"][2, 3] = np.nan

    treated = build(values)
    result = field_of(subtract(treated, treated), "reactivity")
    known = ~np.isnan(field_of(treated, "reactivity"))

    assert known.any(), "nothing was known, so nothing was tested"
    assert np.all(result[known] == 0)
    assert np.isnan(result[~known]).all()


def test_swapping_the_two_inputs_negates_only_the_reactivity(build, subtract, tmp_path):
    treated, untreated = build(everything(seed=4)), build(everything(seed=40))

    forward = subtract(treated, untreated)
    backward = run_subtract(untreated, treated, tmp_path / "backward.h5")

    assert np.array_equal(field_of(backward, "reactivity"),
                          -field_of(forward, "reactivity"))

    for name in set(COMBINED) - {"reactivity"}:
        assert np.array_equal(field_of(backward, name), field_of(forward, name)), name


def test_a_background_of_zeros_leaves_the_treated_run_unchanged(build, subtract):
    treated = build(everything(seed=5))
    nothing = build({"reactivity": 0.0, "error": 0.0})

    output = subtract(treated, nothing)

    for name in ("reactivity", "error", "coverage", "reads/counted"):
        assert np.array_equal(field_of(output, name), field_of(treated, name)), name


def test_two_runs_agree_byte_for_byte(build, tmp_path):
    treated, untreated = build(everything(seed=6)), build(everything(seed=60))

    first = run_subtract(treated, untreated, tmp_path / "first.h5")
    second = run_subtract(treated, untreated, tmp_path / "second.h5")

    assert first.read_bytes() == second.read_bytes()


@pytest.mark.parametrize("n_refs, cap", [(1, 1), (1, 40), (3, 1), (400, 2)])
def test_each_field_follows_its_rule_at_any_shape(build, subtract, n_refs, cap):
    values = {"reactivity": spread("reactivity", n_refs, cap, seed=11),
              "reads/counted": spread("reads/counted", n_refs, cap, seed=12)}
    other = {"reactivity": spread("reactivity", n_refs, cap, seed=21),
             "reads/counted": spread("reads/counted", n_refs, cap, seed=22)}

    treated = build(values, n_refs=n_refs, cap=cap)
    untreated = build(other, n_refs=n_refs, cap=cap)
    output = subtract(treated, untreated)

    for name in ("reactivity", "reads/counted"):
        assert np.array_equal(
            field_of(output, name), combined(name, treated, untreated),
        ), name


# ---------------------------------------------------------------------------
# Inputs that are refused
# ---------------------------------------------------------------------------


def test_inputs_of_different_reference_counts_are_refused(build, tmp_path):
    attempt = try_subtract(build(n_refs=4), build(n_refs=5), tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "references" in attempt.stderr


def test_inputs_of_different_widths_are_refused(build, tmp_path):
    attempt = try_subtract(build(cap=6), build(cap=7), tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "wide" in attempt.stderr


def test_a_file_holding_no_references_is_refused(tmp_path):
    empty = write_output(tmp_path / "empty.h5", n_refs=0, cap=CAP)

    attempt = try_subtract(empty, empty, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "no references" in attempt.stderr


def test_something_that_is_not_hdf5_is_refused(build, tmp_path):
    notes = tmp_path / "notes.txt"
    notes.write_text("months of irreplaceable notes\n")

    attempt = try_subtract(notes, build(), tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert notes.read_text() == "months of irreplaceable notes\n"


@pytest.mark.parametrize("missing", COMBINED)
def test_an_input_missing_any_dataset_is_refused(build, tmp_path, missing):
    attempt = try_subtract(without(build(), missing), build(), tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert missing.rsplit("/", 1)[-1] in attempt.stderr


def test_an_input_whose_datasets_disagree_is_refused(build, tmp_path):
    broken = rewidened(build(), "error", CAP + 3)

    attempt = try_subtract(broken, build(), tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "error" in attempt.stderr


# ---------------------------------------------------------------------------
# What is already at the output path
# ---------------------------------------------------------------------------


def test_an_existing_output_is_not_replaced_without_overwrite(build, subtract):
    treated, untreated = build(), build()
    output = subtract(treated, untreated)
    before = output.read_bytes()

    attempt = try_subtract(treated, untreated, output)

    assert attempt.returncode != 0
    assert "already holds data" in attempt.stderr
    assert output.read_bytes() == before, "the first result is untouched"


def test_overwrite_replaces_an_existing_output(build, subtract):
    """The two runs are given different inputs, so that the values left at the
    path identify which of them wrote it."""
    output = subtract(build({"reactivity": 0.25}), build({"reactivity": 0.0}))

    treated, untreated = build({"reactivity": 0.75}), build({"reactivity": 0.25})
    run_subtract(treated, untreated, output, overwrite=True)

    assert np.all(field_of(output, "reactivity") == np.float32(0.5))


@pytest.mark.parametrize("wrong", ["shape", "missing dataset", "not hdf5"])
def test_a_run_that_refuses_its_inputs_leaves_the_output_intact(build, subtract,
                                                                tmp_path, wrong):
    output = subtract(build(), build())
    before = output.read_bytes()

    def not_hdf5():
        notes = tmp_path / "notes.txt"
        notes.write_text("notes")
        return notes

    bad = {
        "shape": lambda: build(n_refs=N_REFS + 1),
        "missing dataset": lambda: without(build(), RUN_TOTAL),
        "not hdf5": not_hdf5,
    }[wrong]()

    attempt = try_subtract(bad, build(), output, overwrite=True)

    assert attempt.returncode != 0
    assert output.read_bytes() == before


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------


def test_both_inputs_are_required(build, tmp_path):
    given = subprocess.run(
        [str(CMUTS_SUB), "-o", str(tmp_path / "out.h5"), str(build())],
        capture_output=True, text=True,
    )

    assert given.returncode == 2
    assert "expected" in given.stderr
    assert not (tmp_path / "out.h5").exists()


# ---------------------------------------------------------------------------
# End to end
# ---------------------------------------------------------------------------


def test_it_reads_what_cmuts_writes(data, tmp_path):
    """Asserts nothing about any value, only that the run succeeds and leaves
    a file shaped like its inputs."""
    treated = tmp_path / "treated.h5"
    untreated = tmp_path / "untreated.h5"

    run_cmuts(data, treated)
    run_cmuts(data, untreated, min_mapq=30)

    output = run_subtract(treated, untreated, tmp_path / "difference.h5")

    assert datasets_of(output) == datasets_of(treated)
