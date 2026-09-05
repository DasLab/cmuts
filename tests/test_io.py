"""The contracts shared by the programs that read outputs.

cmuts sub, div and norm all read files in the layout cmuts hmm writes, and
write their result in that same layout. Which inputs they accept, what they
leave at the output path, and the form of what they write are one suite of
contracts, run here against each program. The values each program computes are
covered in its own file, as is cmuts hmm's own output path, since its inputs
are alignments and not outputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import pytest

from combining import DIV_RULES, SUB_RULES, expected
from inputs import (
    CAP,
    LARGE_COUNT,
    N_REFS,
    NOTES,
    missing_in_each_input,
    not_hdf5,
    random_fields,
)
from outputs import (
    ALL_FIELDS,
    COUNTED,
    COVERAGE,
    ERROR_FIELDS,
    FIELDS,
    MISMATCH_ERROR,
    FLOAT_FIELDS,
    NORM,
    MISMATCH_RATE,
    RATE_FIELDS,
    SEQUENCE,
    UNMAPPED,
    add_field,
    attributes_of,
    delete_field,
    field_of,
    layout_of,
    outputs_agree,
    read_summary,
    reference_sequences,
    set_field_width,
    shape,
    write_output,
)
from programs import (
    CMUTS_DIV,
    CMUTS_NORM,
    CMUTS_SUB,
    attempt,
    reported_version,
    run_cmuts,
    run_divide,
    run_normalize,
    run_subtract,
    try_divide,
    try_normalize,
    try_subtract,
)

# The error of a ratio rounds at every step of the division and the root, so it
# is compared to a relative tolerance. Every other rule rounds once and is exact.
TOLERANCE = 1e-6


@dataclass(frozen=True)
class Program:
    """One program reading outputs, and how to drive it.

    run and attempt take the inputs as one list, whatever shape the program's
    own command line gives them.
    """

    name: str
    command: tuple
    arity: int        # how many inputs a run here drives it with
    required: tuple   # the datasets an input must hold
    added: tuple      # the datasets the program adds to its output
    rules: dict       # how each field is formed, None where the file says
    rounds: tuple     # the fields compared to a tolerance and not exactly
    run: Callable
    attempt: Callable

    def __str__(self) -> str:
        return self.name


SUBTRACT = Program(
    name="sub", command=CMUTS_SUB, arity=2,
    required=(COVERAGE,) + RATE_FIELDS, added=(), rules=SUB_RULES, rounds=(),
    run=lambda inputs, output, **options: run_subtract(*inputs, output, **options),
    attempt=lambda inputs, output, **options: try_subtract(*inputs, output, **options),
)

DIVIDE = Program(
    name="div", command=CMUTS_DIV, arity=2,
    required=(COVERAGE,) + RATE_FIELDS, added=(), rules=DIV_RULES,
    rounds=ERROR_FIELDS,
    run=lambda inputs, output, **options: run_divide(*inputs, output, **options),
    attempt=lambda inputs, output, **options: try_divide(*inputs, output, **options),
)

NORMALIZE = Program(
    name="norm", command=CMUTS_NORM, arity=1,
    required=(COVERAGE,) + RATE_FIELDS, added=(NORM,), rules=None, rounds=(),
    run=lambda inputs, output, **options: run_normalize(inputs, [output], **options)[0],
    attempt=lambda inputs, output, **options: try_normalize(inputs, [output], **options),
)

READERS = (SUBTRACT, DIVIDE, NORMALIZE)

# The programs forming their output from two inputs by the rules in
# combining.py.
COMBINERS = (SUBTRACT, DIVIDE)

readers = pytest.mark.parametrize("program", READERS, ids=str)
combiners = pytest.mark.parametrize("program", COMBINERS, ids=str)


def inputs_with(program, build, first) -> list:
    """Returns a full set of inputs, with the one under test first and an
    ordinary input in every other position."""
    return [first] + [build() for _ in range(program.arity - 1)]


def written_layout(program, output) -> dict:
    """Returns the layout of an output, less the datasets the program adds."""
    return {name: value for name, value in layout_of(output).items()
            if name not in program.added}


def agrees(program, name, result, expectation) -> bool:
    """Returns whether a field matches its rule: exactly, unless the program
    rounds the field at more than one step."""
    if name in program.rounds:
        return np.allclose(result, expectation, rtol=TOLERANCE, equal_nan=True)

    return np.array_equal(result, expectation)


# ---------------------------------------------------------------------------
# The layout, anchored to what cmuts hmm writes
# ---------------------------------------------------------------------------


def test_the_layout_written_here_is_the_one_cmuts_hmm_writes(data, falsifiable,
                                                             tmp_path):
    """Checks the description in outputs.py against a real cmuts hmm run.

    Compares names, types and widths, which keeps the rate calculation
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


# ---------------------------------------------------------------------------
# Inputs that are refused, and inputs that are read past
# ---------------------------------------------------------------------------


@readers
def test_something_that_is_not_hdf5_is_refused(program, build, tmp_path):
    notes = not_hdf5(tmp_path)
    output = tmp_path / "out.h5"

    failed = program.attempt(inputs_with(program, build, notes), output)

    assert failed.returncode != 0
    assert notes.read_text() == NOTES
    assert not output.exists()


@readers
def test_a_file_holding_no_references_is_refused(program, build, tmp_path):
    empty = write_output(tmp_path / "empty.h5", n_refs=0, cap=CAP)

    failed = program.attempt(inputs_with(program, build, empty), tmp_path / "out.h5")

    assert failed.returncode != 0


@readers
def test_an_input_missing_a_required_dataset_is_refused(program, build, tmp_path):
    for name in program.required:
        for position in range(program.arity):
            inputs = [build() for _ in range(program.arity)]
            inputs[position] = delete_field(build(), name)

            failed = program.attempt(inputs, tmp_path / "out.h5")

            assert failed.returncode != 0, f"{name} missing from input {position}"


@readers
def test_an_input_missing_a_dataset_that_is_not_required_is_skipped(program, build,
                                                                    tmp_path):
    written = 0

    for name in set(ALL_FIELDS) - set(program.required):
        for position in range(program.arity):
            inputs = [build() for _ in range(program.arity)]
            inputs[position] = delete_field(build(), name)

            output = program.run(inputs, tmp_path / f"out{written}.h5")
            written += 1

            assert name not in layout_of(output), f"{name} missing from input {position}"
            assert MISMATCH_RATE in layout_of(output)


@readers
def test_an_input_whose_datasets_disagree_is_refused(program, build, tmp_path):
    broken = set_field_width(build(), MISMATCH_ERROR, CAP + 3)

    failed = program.attempt(inputs_with(program, build, broken), tmp_path / "out.h5")

    assert failed.returncode != 0


@readers
def test_a_dataset_outside_the_layout_is_not_carried(program, build, tmp_path):
    inputs = inputs_with(program, build, add_field(build(), "unknown", np.arange(N_REFS)))

    output = program.run(inputs, tmp_path / "out.h5")

    assert "unknown" not in layout_of(output)
    assert MISMATCH_RATE in layout_of(output)


@combiners
def test_inputs_that_disagree_in_shape_are_refused(program, build, tmp_path):
    rows = program.attempt([build(n_refs=N_REFS + 1), build()], tmp_path / "rows.h5")
    columns = program.attempt([build(cap=CAP + 1), build()], tmp_path / "columns.h5")

    assert rows.returncode != 0
    assert columns.returncode != 0


@combiners
def test_inputs_that_disagree_on_the_sequence_are_refused(program, build, tmp_path):
    changed = reference_sequences(N_REFS, CAP)
    changed[0, 0] = (changed[0, 0] + 1) % 4

    failed = program.attempt([build(), build({SEQUENCE: changed})], tmp_path / "out.h5")

    assert failed.returncode != 0


# ---------------------------------------------------------------------------
# An existing file at the output path
# ---------------------------------------------------------------------------


@readers
def test_an_existing_output_is_not_replaced_without_overwrite(program, build, tmp_path):
    inputs = [build() for _ in range(program.arity)]
    output = program.run(inputs, tmp_path / "out.h5")
    before = output.read_bytes()

    failed = program.attempt(inputs, output)

    assert failed.returncode != 0
    assert output.read_bytes() == before


@readers
def test_a_file_that_is_not_an_output_is_not_replaced_without_overwrite(program, build,
                                                                        tmp_path):
    notes = not_hdf5(tmp_path)

    failed = program.attempt([build() for _ in range(program.arity)], notes)

    assert failed.returncode != 0
    assert notes.read_text() == NOTES


@readers
def test_an_empty_file_is_replaced_without_overwrite(program, build, tmp_path):
    output = tmp_path / "reserved.h5"
    output.touch()

    program.run([build() for _ in range(program.arity)], output)

    assert MISMATCH_RATE in layout_of(output)


@readers
def test_overwrite_replaces_an_existing_output(program, build, tmp_path):
    """The two runs are given different inputs, so that the values left at the
    path identify which of them wrote it."""
    output = program.run([build(random_fields(seed=70 + i))
                          for i in range(program.arity)], tmp_path / "out.h5")
    before = output.read_bytes()

    second = [build(random_fields(seed=80 + i)) for i in range(program.arity)]
    program.run(second, output, overwrite=True)
    separate = program.run(second, tmp_path / "separate.h5")

    assert output.read_bytes() != before
    assert outputs_agree(output, separate)


@readers
def test_a_run_that_refuses_its_inputs_leaves_the_output_intact(program, build,
                                                                tmp_path):
    output = program.run([build() for _ in range(program.arity)], tmp_path / "out.h5")
    before = output.read_bytes()

    refused = {
        "missing dataset": delete_field(build(), MISMATCH_RATE),
        "not hdf5": not_hdf5(tmp_path),
    }

    for wrong, bad in refused.items():
        failed = program.attempt(inputs_with(program, build, bad), output,
                                 overwrite=True)

        assert failed.returncode != 0, wrong
        assert output.read_bytes() == before, wrong


# ---------------------------------------------------------------------------
# The written form
# ---------------------------------------------------------------------------


@readers
def test_the_output_is_shaped_and_typed_like_its_input(program, build, tmp_path):
    inputs = [build(random_fields(seed=95 + i)) for i in range(program.arity)]

    output = program.run(inputs, tmp_path / "out.h5")

    assert set(layout_of(output)) == set(layout_of(inputs[0])) | set(program.added)
    assert written_layout(program, output) == layout_of(inputs[0])


@readers
def test_the_result_names_the_program_that_wrote_it(program, build, tmp_path):
    output = program.run([build() for _ in range(program.arity)], tmp_path / "out.h5")

    assert attributes_of(output)["program"] == " ".join(program.command)
    assert attributes_of(output)["version"] == reported_version(program.command)


@readers
def test_two_runs_agree_byte_for_byte(program, build, tmp_path):
    inputs = [build(random_fields(seed=90 + i)) for i in range(program.arity)]

    first = program.run(inputs, tmp_path / "first.h5")
    second = program.run(inputs, tmp_path / "second.h5")

    assert first.read_bytes() == second.read_bytes()


@readers
def test_the_columns_past_a_reference_stay_nan(program, build, tmp_path):
    lengths = [6, 4, 2, 1]
    rows = np.full((N_REFS, CAP), np.float32(0.25))

    for row, length in enumerate(lengths):
        rows[row, length:] = np.nan

    inputs = [build({MISMATCH_RATE: rows}) for _ in range(program.arity)]
    output = program.run(inputs, tmp_path / "out.h5")
    result = field_of(output, MISMATCH_RATE)

    for row, length in enumerate(lengths):
        assert not np.isnan(result[row, :length]).any(), row
        assert np.isnan(result[row, length:]).all(), row


# ---------------------------------------------------------------------------
# The rules the combiners share in form
# ---------------------------------------------------------------------------


@combiners
@pytest.mark.parametrize("name", ALL_FIELDS)
def test_each_field_follows_its_rule(program, build, tmp_path, name):
    inputs = [build(random_fields(seed=1), unmapped=17),
              build(random_fields(seed=50), unmapped=3)]

    output = program.run(inputs, tmp_path / "out.h5")

    assert agrees(program, name, field_of(output, name),
                  expected(program.rules, name, *inputs)), name


@combiners
@pytest.mark.parametrize("n_refs, cap", [(1, 1), (1, 40), (3, 1), (400, 2)])
def test_each_field_follows_its_rule_at_any_shape(program, build, tmp_path,
                                                  n_refs, cap):
    shaped = dict(n_refs=n_refs, cap=cap)
    inputs = [build(random_fields(seed=11, **shaped), unmapped=17, **shaped),
              build(random_fields(seed=21, **shaped), unmapped=3, **shaped)]

    output = program.run(inputs, tmp_path / "out.h5")

    for name in ALL_FIELDS:
        assert agrees(program, name, field_of(output, name),
                      expected(program.rules, name, *inputs)), name


@combiners
def test_counts_stay_whole_and_exact(program, build, tmp_path):
    inputs = [build({COUNTED: LARGE_COUNT}, unmapped=LARGE_COUNT),
              build({COUNTED: 1}, unmapped=1)]

    output = program.run(inputs, tmp_path / "out.h5")

    assert np.all(field_of(output, COUNTED) == LARGE_COUNT + 1)
    assert field_of(output, UNMAPPED) == LARGE_COUNT + 1
    assert field_of(output, COUNTED).dtype == np.dtype("u8")


@combiners
@pytest.mark.parametrize("name", FLOAT_FIELDS)
def test_a_value_either_input_lacks_is_missing_from_the_output(program, build,
                                                               tmp_path, name):
    left, right = missing_in_each_input(N_REFS, CAP)

    # The error of a ratio uses its channel's rates, so every rate is set to a
    # known value and the field under test is the only one with a value missing.
    known = dict.fromkeys(RATE_FIELDS, 0.5)

    output = program.run([build(known | {name: left}), build(known | {name: right})],
                         tmp_path / "out.h5")
    result = field_of(output, name)

    assert np.array_equal(np.isnan(result), np.isnan(left) | np.isnan(right))
    assert not np.isnan(result[0]).any()


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------


@combiners
def test_every_input_is_required(program, build, tmp_path):
    given = attempt([*program.command, "-o", tmp_path / "out.h5", build()])

    assert given.returncode == 2
    assert not (tmp_path / "out.h5").exists()


# ---------------------------------------------------------------------------
# End to end
# ---------------------------------------------------------------------------


@readers
def test_the_program_reads_what_cmuts_hmm_writes(program, data, falsifiable, tmp_path):
    """Asserts that the run succeeds and leaves a file shaped like its
    inputs."""
    counted = tmp_path / "counted.h5"
    summary = read_summary(run_cmuts(data, counted))

    falsifiable(summary.rows > 0)

    inputs = [counted]
    if program.arity == 2:
        inputs.append(run_cmuts(data, tmp_path / "stricter.h5", min_mapq=30))

    output = program.run(inputs, tmp_path / "combined.h5")

    assert written_layout(program, output) == layout_of(counted)
