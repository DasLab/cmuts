"""What a value of NaN means in the output.

The arrays with a row per reference are rectangular and the references are not,
so every row runs to the width the longest reference needs. NaN marks the two
cases this leaves: a column past what its own reference needs, and a reference
that no read was aligned to.

A zero means something else: the position was reached and nothing was counted
there. Anything summing these arrays depends on the distinction, so the meaning
of each value is fixed here and not left to the writer.
"""

import itertools
from dataclasses import dataclass

import numpy as np
import pytest

from oracle import reference_lengths, references_with_reads, rows_by_name
from outputs import (
    COLUMN_COUNT_FIELDS,
    COVERAGE,
    ERROR,
    PADDED_FIELDS,
    BY_NAME,
    PER_BASE_FIELDS,
    PER_REFERENCE_FIELDS,
    REACTIVITY,
    REJECTED,
    extent,
    field_of,
    fields_of,
)
from programs import run_cmuts

# Higher than every mapping quality the generator writes. A read of
# unavailable quality is rejected at any threshold, so every read is rejected.
REJECTS_EVERYTHING = 61

# Far enough apart that the middle one is not the only step, and starting below
# one whole observation so that no position is discarded to begin with.
DEPTHS = (0, 1, 5)


@dataclass(frozen=True)
class ScoredRun:
    """A run, and what a test needs in order to read its rows: where each
    reference was written, how far it reaches, and whether any read aligned to
    it."""

    output: object
    lengths: dict
    row_of: dict
    reached: set

    @property
    def widest(self) -> int:
        return max(self.lengths.values())

    @property
    def shorter(self) -> list:
        """The references that the widest one leaves padding after."""
        return [name for name, length in self.lengths.items() if length < self.widest]

    @property
    def missing(self) -> list:
        """The references that no read aligned to."""
        return [name for name in self.row_of if name not in self.reached]


@pytest.fixture
def scored_run(tmp_path):
    """Returns a function that runs one dataset and gathers what its output is
    read against. Each call writes a file of its own, so one test may make
    several."""
    written = itertools.count()

    def run(data, **options):
        output = tmp_path / f"{data.name}{next(written)}.h5"
        run_cmuts(data, output, **options)

        return ScoredRun(
            output=output,
            lengths=reference_lengths(data.fasta),
            row_of=rows_by_name(data.fasta),
            reached=references_with_reads(data.bam),
        )

    return run


# ---------------------------------------------------------------------------
# A reference shorter than the row it is written to
# ---------------------------------------------------------------------------


def test_positions_past_a_reference_are_nan(data, scored_run, falsifiable):
    run = scored_run(data)

    falsifiable(len(run.shorter) > 0)

    fields = fields_of(run.output, PER_BASE_FIELDS)

    for reference in run.shorter:
        for field, values in fields.items():
            tail = values[run.row_of[reference]][run.lengths[reference]:]
            pad = np.full(tail.shape, BY_NAME[field].pad, dtype=tail.dtype)

            assert np.array_equal(tail, pad, equal_nan=tail.dtype.kind == "f"), \
                f"{field}: {reference} is padded with {tail[:4]}"


def test_positions_within_a_reference_are_never_nan(data, scored_run, falsifiable):
    run = scored_run(data)

    falsifiable(len(run.reached) > 0)

    fields = fields_of(run.output, PADDED_FIELDS)

    for reference in run.reached:
        for field, values in fields.items():
            within = values[run.row_of[reference]][:run.lengths[reference]]

            assert not np.isnan(within).any(), f"{field}: {reference} has a hole in it"


# ---------------------------------------------------------------------------
# A reference that no read was aligned to
# ---------------------------------------------------------------------------


def test_a_reference_with_no_reads_is_zero_over_its_own_bases(data, scored_run,
                                                              falsifiable):
    run = scored_run(data)

    falsifiable(len(run.missing) > 0)

    counted = fields_of(run.output, COLUMN_COUNT_FIELDS)
    holes = fields_of(run.output, PADDED_FIELDS)
    totals = fields_of(run.output, PER_REFERENCE_FIELDS)

    for reference in run.missing:
        row = run.row_of[reference]
        length = run.lengths[reference]

        for field, values in counted.items():
            columns = extent(field, length, values.shape[1])

            assert (values[row][:columns] == 0).all(), \
                f"{field}: {reference} has no reads and is not zero"

        for field, values in holes.items():
            assert np.isnan(values[row][length:]).all(), \
                f"{field}: {reference} has padding that is not NaN"

        # A count of reads is zero where no read arrived, so no NaN appears
        # here.
        for field, values in totals.items():
            assert values[row] == 0, \
                f"{field}: {reference} has no reads but is not zero"


def test_an_uncovered_reference_of_full_length_holds_no_nan(data, scored_run,
                                                            falsifiable):
    run = scored_run(data)
    full_length = [name for name in run.missing if run.lengths[name] == run.widest]

    falsifiable(len(full_length) > 0)

    holes = fields_of(run.output, PADDED_FIELDS)
    counted = fields_of(run.output, COLUMN_COUNT_FIELDS)

    for reference in full_length:
        row = run.row_of[reference]

        for field, values in holes.items():
            assert not np.isnan(values[row]).any(), f"{field}: {reference} holds a NaN"

        for field, values in counted.items():
            assert (values[row] == 0).all(), f"{field}: {reference} is not zero"


def test_a_reference_whose_reads_were_all_rejected_is_zero(data, scored_run,
                                                           falsifiable):
    run = scored_run(data, min_mapq=REJECTS_EVERYTHING)

    falsifiable(len(run.reached) > 0)

    holes = fields_of(run.output, PADDED_FIELDS)
    counted = fields_of(run.output, COLUMN_COUNT_FIELDS)
    totals = fields_of(run.output, PER_REFERENCE_FIELDS)

    for reference in run.reached:
        row = run.row_of[reference]
        length = run.lengths[reference]

        for field, values in holes.items():
            assert not np.isnan(values[row][:length]).any(), \
                f"{field}: {reference} went to NaN"

        for field, values in counted.items():
            columns = extent(field, length, values.shape[1])

            assert (values[row][:columns] == 0).all(), \
                f"{field}: {reference} counted something"

        for field, values in totals.items():
            assert values[row] == 0 or field == REJECTED, \
                f"{field}: {reference} is not zero"


# ---------------------------------------------------------------------------
# What NaN means in a rate
# ---------------------------------------------------------------------------


def test_a_rate_is_known_wherever_its_error_is(data, scored_run, falsifiable):
    """Separates the two reasons a rate is NaN by the coverage, which is NaN
    outside a reference and nowhere else."""
    run = scored_run(data)

    fields = fields_of(run.output, (REACTIVITY, ERROR, COVERAGE))
    reactivity, error, coverage = fields[REACTIVITY], fields[ERROR], fields[COVERAGE]

    finite = reactivity[~np.isnan(reactivity)]

    falsifiable(finite.size > 0)

    assert np.array_equal(np.isnan(reactivity), np.isnan(error)), \
        "the rate and its error are not NaN at the same positions"

    # Padding is not part of any reference, so no rate is had there either.
    assert np.isnan(reactivity[np.isnan(coverage)]).all(), \
        "a rate outside a reference is not NaN"

    assert ((finite >= 0) & (finite <= 1)).all(), \
        f"a rate of {finite.min()} to {finite.max()} lies outside zero to one"


def test_raising_the_minimum_depth_only_discards_rates(data, scored_run, falsifiable):
    missing = {}

    for depth in DEPTHS:
        run = scored_run(data, min_depth=depth)
        missing[depth] = np.isnan(field_of(run.output, REACTIVITY))

    falsifiable(missing[DEPTHS[-1]].sum() > missing[DEPTHS[0]].sum())

    for lower, higher in itertools.pairwise(DEPTHS):
        assert (missing[higher] >= missing[lower]).all(), \
            f"depth {higher} recovered a rate that depth {lower} had not"


def test_an_error_at_a_full_read_of_depth_is_at_most_a_half(data, scored_run,
                                                            falsifiable):
    """Runs at --min-depth 1, below which the standard error of a proportion
    can exceed a half."""
    run = scored_run(data, min_depth=1)

    error = field_of(run.output, ERROR)
    finite = error[~np.isnan(error)]

    falsifiable(finite.size > 0)

    assert (finite <= 0.5).all(), f"an error of {finite.max()} over a whole read"
