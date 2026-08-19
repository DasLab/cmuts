"""Measuring reactivity rates against a known structure.

A row of the table depends on the values in the input file, the sequence and
the pairing, and on nothing else, so all three are written by hand. scoring.py
builds them and computes what each row should hold.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from inputs import CAP, N_REFS, not_hdf5, random_values
from outputs import COVERAGE, REACTIVITY, read_summary
from programs import CMUTS_SCORE, attempt, run_cmuts, run_score, try_score
from scoring import (
    COLUMNS,
    alternating,
    expected,
    kept,
    order_of,
    rows_of,
    sequences_of,
    write_fasta,
    write_structures,
)

# auroc and auprc are written to five decimals and a mean to six significant
# digits, so a comparison against an oracle in double allows for the rounding.
TOLERANCE = 1e-5

NAME = "ref0"
SEQUENCE = "ACGUAC"


@pytest.fixture
def score(tmp_path):
    """Returns a function that scores a file against a structure written for
    it. Each pair of files is named separately, so one test may score
    several."""
    written = itertools.count()

    def run(path, sequences, pairings, **options):
        which = next(written)
        fasta = write_fasta(tmp_path / f"refs{which}.fasta", sequences)
        structures = write_structures(tmp_path / f"refs{which}.db", pairings)

        return run_score(path, fasta, structures, **options)

    return run


def one_row(values):
    """Puts the values in the first row of a file, leaving every other row
    unmeasured."""
    rows = np.full((N_REFS, CAP), np.nan, dtype=np.float32)
    rows[0] = values

    return rows


def only(table: dict) -> dict:
    """The single row a table of one reference holds."""
    (row,) = table.values()

    return row


# ---------------------------------------------------------------------------
# The table
# ---------------------------------------------------------------------------


def test_the_header_names_every_column(build, score):
    written = score(build({REACTIVITY: one_row([0, 1, 0, 1, 0, 1])}),
                    {NAME: SEQUENCE}, {NAME: ".(.(.("})

    assert written.splitlines()[0] == ",".join(COLUMNS)


def test_the_rows_follow_the_fasta(build, score):
    names = ["a", "b", "c"]
    values = np.tile([0, 1, 0, 1, 0, 1], (N_REFS, 1)).astype(np.float32)

    written = score(build({REACTIVITY: values}),
                    {name: SEQUENCE for name in names},
                    {name: ".(.(.(" for name in reversed(names)})

    assert order_of(written) == names


def test_a_reference_with_no_structure_is_left_out(build, score):
    values = np.tile([0, 1, 0, 1, 0, 1], (N_REFS, 1)).astype(np.float32)

    written = score(build({REACTIVITY: values}),
                    {"a": SEQUENCE, "b": SEQUENCE}, {"b": ".(.(.("})

    assert order_of(written) == ["b"]


# ---------------------------------------------------------------------------
# The metrics
# ---------------------------------------------------------------------------


def test_a_perfect_ranking_scores_one(build, score):
    """The unpaired bases carry the three highest values."""
    table = rows_of(score(build({REACTIVITY: one_row([9, 1, 8, 2, 7, 3])}),
                          {NAME: SEQUENCE}, {NAME: ".(.(.("}))

    assert only(table)["auroc"] == 1.0
    assert only(table)["auprc"] == 1.0


def test_a_reversed_ranking_scores_zero(build, score):
    table = rows_of(score(build({REACTIVITY: one_row([1, 9, 2, 8, 3, 7])}),
                          {NAME: SEQUENCE}, {NAME: ".(.(.("}))

    assert only(table)["auroc"] == 0.0


def test_one_value_everywhere_scores_a_half(build, score):
    """Every pair is a tie, and a tie counts as half."""
    table = rows_of(score(build({REACTIVITY: one_row([0.4] * CAP)}),
                          {NAME: SEQUENCE}, {NAME: ".(.(.("}))

    assert only(table)["auroc"] == 0.5


@pytest.mark.parametrize("rounding", ["distinct", "tied"])
def test_every_metric_matches_the_oracle(build, score, rounding):
    """Ties are scored by a rule of their own, so the values are run both as
    they come and rounded until they repeat."""
    values = random_values(REACTIVITY, seed=51)

    if rounding == "tied":
        values = np.round(values)

    pairing = alternating(CAP)
    path = build({REACTIVITY: values})
    table = rows_of(score(path, {NAME: SEQUENCE}, {NAME: pairing}))

    wanted = expected(*kept(path, 0, SEQUENCE, pairing))

    for column, value in wanted.items():
        assert np.isclose(only(table)[column], value, rtol=TOLERANCE), column


def test_the_means_are_of_each_class(build, score):
    table = rows_of(score(build({REACTIVITY: one_row([1, 3, 1, 3, 1, 3])}),
                          {NAME: SEQUENCE}, {NAME: ".(.(.("}))

    assert only(table)["mean_unpaired"] == 1.0
    assert only(table)["mean_paired"] == 3.0


def test_a_rise_that_keeps_the_order_leaves_the_auroc_alone(build, score):
    """auroc reads the ranking and not the values, so scaling every rate leaves
    it where it was."""
    values = random_values(REACTIVITY, seed=52)
    pairing = alternating(CAP)

    first = rows_of(score(build({REACTIVITY: values}), {NAME: SEQUENCE},
                          {NAME: pairing}))
    second = rows_of(score(build({REACTIVITY: values * 3 + 1}), {NAME: SEQUENCE},
                           {NAME: pairing}))

    assert np.isclose(only(first)["auroc"], only(second)["auroc"], rtol=TOLERANCE)


def test_swapping_the_classes_complements_the_auroc(build, score):
    values = random_values(REACTIVITY, seed=53)
    pairing = alternating(CAP)
    swapped = "".join("(" if mark == "." else "." for mark in pairing)

    first = rows_of(score(build({REACTIVITY: values}), {NAME: SEQUENCE},
                          {NAME: pairing}))
    second = rows_of(score(build({REACTIVITY: values}), {NAME: SEQUENCE},
                           {NAME: swapped}))

    assert np.isclose(only(first)["auroc"], 1 - only(second)["auroc"],
                      rtol=TOLERANCE)


# ---------------------------------------------------------------------------
# What is scored
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("kind", ["()", "[]", "{}", "<>"])
def test_every_bracket_counts_as_paired(build, score, kind):
    pairing = "." + kind[0] + "." + kind[1] + ".."
    table = rows_of(score(build({REACTIVITY: one_row([1, 2, 3, 4, 5, 6])}),
                          {NAME: SEQUENCE}, {NAME: pairing}))

    assert only(table)["paired"] == 2
    assert only(table)["unpaired"] == 4


def test_a_mark_that_is_neither_dot_nor_bracket_is_not_scored(build, score):
    table = rows_of(score(build({REACTIVITY: one_row([1, 2, 3, 4, 5, 6])}),
                          {NAME: SEQUENCE}, {NAME: ".((--)"}))

    assert only(table)["paired"] + only(table)["unpaired"] == 4


def test_a_base_outside_the_bases_is_not_scored(build, score):
    """The sequence is ACGUAC, so asking for A and C leaves four positions."""
    table = rows_of(score(build({REACTIVITY: one_row([1, 2, 3, 4, 5, 6])}),
                          {NAME: SEQUENCE}, {NAME: ".(.(.("}, bases="A,C"))

    assert only(table)["paired"] + only(table)["unpaired"] == 4


def test_a_reference_written_as_dna_is_scored_at_u(build, score):
    """T and U name one base, so a DNA reference is scored at its T."""
    table = rows_of(score(build({REACTIVITY: one_row([1, 2, 3, 4, 5, 6])}),
                          {NAME: "ACGTAT"}, {NAME: ".(.((."}, bases="U"))

    assert only(table)["paired"] == 1
    assert only(table)["unpaired"] == 1


def test_a_position_below_the_coverage_floor_is_not_scored(build, score):
    path = build({REACTIVITY: one_row([1, 2, 3, 4, 5, 6]),
                  COVERAGE: one_row([100, 100, 100, 1, 1, 1])})

    table = rows_of(score(path, {NAME: SEQUENCE}, {NAME: ".(.(.("},
                          min_coverage=50))

    assert only(table)["paired"] + only(table)["unpaired"] == 3


def test_a_position_with_no_reactivity_is_not_scored(build, score):
    table = rows_of(score(build({REACTIVITY: one_row([1, 2, np.nan, 4, 5, 6])}),
                          {NAME: SEQUENCE}, {NAME: ".(.(.("}))

    assert only(table)["paired"] + only(table)["unpaired"] == 5


def test_a_record_may_carry_its_sequence(build, tmp_path):
    """A dot bracket file often holds the sequence before the pairing, which is
    passed over."""
    path = build({REACTIVITY: one_row([9, 1, 8, 2, 7, 3])})
    fasta = write_fasta(tmp_path / "refs.fasta", {NAME: SEQUENCE})
    structures = write_structures(tmp_path / "refs.db", {NAME: ".(.(.("},
                                  sequences={NAME: SEQUENCE})

    table = rows_of(run_score(path, fasta, structures))

    assert only(table)["auroc"] == 1.0


# ---------------------------------------------------------------------------
# What is refused
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("pairing", ["......", "(((((("])
def test_a_reference_of_one_class_is_left_out(build, tmp_path, pairing):
    """A ranking needs both classes, so a reference holding one is not scored,
    and a run left with no other reference fails."""
    path = build({REACTIVITY: one_row([1, 2, 3, 4, 5, 6])})
    fasta = write_fasta(tmp_path / "refs.fasta", {NAME: SEQUENCE})
    structures = write_structures(tmp_path / "refs.db", {NAME: pairing})

    failed = try_score(path, fasta, structures)

    assert failed.returncode == 1
    assert failed.stdout == ""


def test_a_structure_of_another_length_is_not_scored(build, tmp_path):
    path = build({REACTIVITY: one_row([1, 2, 3, 4, 5, 6])})
    fasta = write_fasta(tmp_path / "refs.fasta", {NAME: SEQUENCE})
    structures = write_structures(tmp_path / "refs.db", {NAME: ".(.("})

    failed = try_score(path, fasta, structures)

    assert failed.returncode == 1
    assert NAME in failed.stderr


def test_a_structure_longer_than_a_line_is_refused(build, tmp_path):
    """A record longer than the reader's line is named as such, and not
    reported as a pairing with no name."""
    length = 8192
    path = build({REACTIVITY: np.zeros((N_REFS, length), dtype=np.float32)},
                 cap=length)
    fasta = write_fasta(tmp_path / "refs.fasta", {NAME: "A" * length})
    structures = write_structures(tmp_path / "refs.db",
                                  {NAME: alternating(length)})

    failed = try_score(path, fasta, structures)

    assert failed.returncode == 1
    assert "too long" in failed.stderr


def test_a_fasta_longer_than_the_file_is_refused(build, tmp_path):
    path = build({REACTIVITY: np.zeros((1, CAP), dtype=np.float32)}, n_refs=1)
    fasta = write_fasta(tmp_path / "refs.fasta", {"a": SEQUENCE, "b": SEQUENCE})
    structures = write_structures(tmp_path / "refs.db",
                                  {"a": ".(.(.(", "b": ".(.(.("})

    assert try_score(path, fasta, structures).returncode == 1


def test_an_input_that_is_not_an_output_is_refused(tmp_path):
    fasta = write_fasta(tmp_path / "refs.fasta", {NAME: SEQUENCE})
    structures = write_structures(tmp_path / "refs.db", {NAME: ".(.(.("})

    assert try_score(not_hdf5(tmp_path), fasta, structures).returncode == 1


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("missing", ["fasta", "structures", "input"])
def test_every_required_argument_is_required(build, tmp_path, missing):
    fasta = write_fasta(tmp_path / "refs.fasta", {NAME: SEQUENCE})
    structures = write_structures(tmp_path / "refs.db", {NAME: ".(.(.("})
    given = {"fasta": ["-f", fasta], "structures": ["-s", structures],
             "input": [build()]}

    words = [word for name, part in given.items() if name != missing
             for word in part]

    assert attempt([*CMUTS_SCORE, *words]).returncode == 2


# ---------------------------------------------------------------------------
# End to end
# ---------------------------------------------------------------------------


def test_cmuts_score_reads_what_cmuts_hmm_writes(data, falsifiable, tmp_path):
    """A reference with no finite reactivity is scored by nothing, and a
    dataset whose references are all such leaves nothing to score, so the run
    is allowed to refuse it."""
    rates = run_cmuts(data, tmp_path / "rates.h5")
    sequences = sequences_of(data.fasta)
    structures = write_structures(
        tmp_path / "refs.db",
        {name: alternating(len(sequence)) for name, sequence in sequences.items()})

    finished = try_score(rates, data.fasta, structures)

    falsifiable(finished.returncode == 0 and read_summary(rates).rows > 0)

    if finished.returncode != 0:
        return

    table = rows_of(finished.stdout)

    assert set(table) <= set(sequences)
    assert all(row["paired"] > 0 and row["unpaired"] > 0 for row in table.values())
