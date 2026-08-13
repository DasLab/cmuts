"""What a value of NaN means in the output.

The arrays with a row per reference are rectangular and the references are not,
so every row runs to the width the longest reference needs. NaN is what tells
the difference: a column past what its own reference needs, and a reference no
read was ever aligned to.

A zero means something else: the position was reached and nothing was counted
there. Anything summing these arrays depends on the distinction, so what each
value means is fixed here rather than left to the writer.
"""

import itertools
from dataclasses import dataclass

import h5py
import numpy as np
import pytest

from datasets import DATASETS
from outputs import BY_LENGTH, RATES
from support import (
    arrays_of, reference_lengths, references_with_reads, rows_by_name, run_cmuts,
)

# Higher than any mapping quality a read can carry, so every read is rejected.
REJECTS_EVERYTHING = 61

# Enough apart that the middle one is not the only step, and starting below the
# first whole observation so that nothing is discarded at all to begin with.
DEPTHS = (0, 1, 5)


def rectangular(output):
    """The arrays with a row per reference, indexed by position or by length."""
    return {k: d[:] for k, d in arrays_of(output).items() if d.ndim == 2}


def counts(output):
    """The rows holding counts, whether a reference's own length bounds them or
    a read length does. These are the ones a zero is meaningful in."""
    return {k: v for k, v in rectangular(output).items() if k not in RATES}


def per_base(output):
    """The arrays a reference's own length bounds, which are the padded ones."""
    return {k: v for k, v in rectangular(output).items() if k not in BY_LENGTH}


def padded(output):
    """The arrays a NaN is meaningful in: bounded by a reference and holding
    counts.

    A row indexed by read length runs to its full width and holds whole
    numbers, so no NaN can arise there. The rates carry one at any position
    failing --min-depth as well, so NaN in them does not mean padding on its
    own.
    """
    return {k: v for k, v in per_base(output).items() if k not in RATES}


def row_extent(field, ref_len, width):
    """Columns of a row that hold data; the rest is padding."""
    return width if field in BY_LENGTH else ref_len


def per_reference(output):
    """The arrays with one value for each reference."""
    return {k: d[:] for k, d in arrays_of(output).items() if d.ndim == 1}


@dataclass(frozen=True)
class Scored:
    """A run, and what a test needs to read its rows: where each reference was
    written, how far it reaches, and whether any read arrived on it."""

    data: object
    output: object
    lengths: dict
    row_of: dict
    reached: set

    @property
    def widest(self) -> int:
        return max(self.lengths.values())

    @property
    def shorter(self) -> list:
        """References the widest one leaves padding after."""
        return [name for name, length in self.lengths.items() if length < self.widest]

    @property
    def missing(self) -> list:
        """References no read reached."""
        return [name for name in self.row_of if name not in self.reached]


@pytest.fixture
def scored(datasets, tmp_path):
    """Runs one dataset, gathering what its output is read against. Each run
    writes a file of its own, so one test may make several."""
    written = itertools.count()

    def run(name, **options):
        data = datasets(name)
        output = tmp_path / f"{name}{next(written)}.h5"
        run_cmuts(data, output, **options)

        return Scored(
            data=data,
            output=output,
            lengths=reference_lengths(data.fasta),
            row_of=rows_by_name(data.fasta),
            reached=references_with_reads(data.bam),
        )

    return run


# ---------------------------------------------------------------------------
# A reference shorter than the row it is written to
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_positions_past_a_reference_are_nan(scored, checked, name):
    run = scored(name)

    with h5py.File(run.output, "r") as handle:
        fields = per_base(handle)

    for reference in checked(run.shorter):
        for field, values in fields.items():
            tail = values[run.row_of[reference]][run.lengths[reference]:]

            assert np.isnan(tail).all(), \
                f"{field}: {reference} is padded with {tail[:4]}"


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_positions_within_a_reference_are_never_nan(scored, checked, name):
    run = scored(name)

    with h5py.File(run.output, "r") as handle:
        fields = padded(handle)

    for reference in checked(run.reached):
        for field, values in fields.items():
            within = values[run.row_of[reference]][:run.lengths[reference]]

            assert not np.isnan(within).any(), f"{field}: {reference} has a hole in it"


# ---------------------------------------------------------------------------
# A reference no read was aligned to
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_reference_with_no_reads_is_zero_over_its_own_bases(scored, checked, name):
    run = scored(name)

    with h5py.File(run.output, "r") as handle:
        counted, holes, totals = counts(handle), padded(handle), per_reference(handle)

    for reference in checked(run.missing):
        row = run.row_of[reference]
        length = run.lengths[reference]

        for field, values in counted.items():
            extent = row_extent(field, length, values.shape[1])

            assert (values[row][:extent] == 0).all(), \
                f"{field}: {reference} has no reads and is not zero"

        for field, values in holes.items():
            assert np.isnan(values[row][length:]).all(), \
                f"{field}: {reference} has padding that is not NaN"

        # A count of reads is zero where no read arrived: NaN would say nothing
        # that zero does not.
        for field, values in totals.items():
            assert values[row] == 0, \
                f"{field}: {reference} has no reads but is not zero"


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_an_uncovered_reference_of_full_length_holds_no_nan(scored, checked, name):
    run = scored(name)
    full_length = [name for name in run.missing if run.lengths[name] == run.widest]

    with h5py.File(run.output, "r") as handle:
        holes, counted = padded(handle), counts(handle)

    for reference in checked(full_length):
        row = run.row_of[reference]

        for field, values in holes.items():
            assert not np.isnan(values[row]).any(), f"{field}: {reference} holds a NaN"

        for field, values in counted.items():
            assert (values[row] == 0).all(), f"{field}: {reference} is not zero"


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_reference_whose_reads_were_all_rejected_is_zero(scored, checked, name):
    run = scored(name, min_mapq=REJECTS_EVERYTHING)

    with h5py.File(run.output, "r") as handle:
        holes, counted, totals = padded(handle), counts(handle), per_reference(handle)

    for reference in checked(run.reached):
        row = run.row_of[reference]
        length = run.lengths[reference]

        for field, values in holes.items():
            assert not np.isnan(values[row][:length]).any(), \
                f"{field}: {reference} went to NaN"

        for field, values in counted.items():
            extent = row_extent(field, length, values.shape[1])

            assert (values[row][:extent] == 0).all(), \
                f"{field}: {reference} counted something"

        for field, values in totals.items():
            assert values[row] == 0 or field == "reads/rejected", \
                f"{field}: {reference} is not zero"


# ---------------------------------------------------------------------------
# What NaN means in a rate
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_rate_is_known_wherever_its_error_is(scored, checked, name):
    """Separates the two reasons a rate is NaN by the coverage, which is NaN
    outside a reference only."""
    run = scored(name)

    with h5py.File(run.output, "r") as handle:
        reactivity = handle["reactivity"][:]
        error = handle["error"][:]
        coverage = handle["coverage"][:]

    assert np.array_equal(np.isnan(reactivity), np.isnan(error)), \
        "the rate and its error disagree about what is known"

    # Padding is not part of any reference, so no rate is had there either.
    assert np.isnan(reactivity[np.isnan(coverage)]).all(), \
        "a rate outside a reference is not NaN"

    finite = checked(reactivity[~np.isnan(reactivity)])

    assert ((finite >= 0) & (finite <= 1)).all(), \
        f"a rate of {finite.min()} to {finite.max()} lies outside zero to one"


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_raising_the_minimum_depth_only_discards_rates(scored, checked, name):
    missing = {}

    for depth in DEPTHS:
        run = scored(name, min_depth=depth)

        with h5py.File(run.output, "r") as handle:
            missing[depth] = np.isnan(handle["reactivity"][:])

    for lower, higher in itertools.pairwise(DEPTHS):
        assert (missing[higher] >= missing[lower]).all(), \
            f"depth {higher} recovered a rate that depth {lower} had not"

    checked(missing[DEPTHS[-1]].sum() - missing[DEPTHS[0]].sum())


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_an_error_at_a_full_read_of_depth_is_at_most_a_half(scored, checked, name):
    """Runs at --min-depth 1, below which the standard error of a proportion
    can exceed a half."""
    run = scored(name, min_depth=1)

    with h5py.File(run.output, "r") as handle:
        error = handle["error"][:]

    finite = checked(error[~np.isnan(error)])

    assert (finite <= 0.5).all(), f"an error of {finite.max()} over a whole read"
