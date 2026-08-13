"""What a value of NaN means in the output.

The arrays with a row per reference are rectangular and the references are not,
so every row runs to the width the longest reference needs. NaN is what tells
the difference: a column past what its own reference needs, and a reference no
read was ever aligned to.

A zero means something else: the position was reached and nothing was counted
there. Anything summing these arrays depends on the distinction, so what each
value means is fixed here rather than left to the writer.
"""

import h5py
import numpy as np
import pytest

from support import (
    RATES, datasets_of, references_with_reads, rows_by_name, run_cmuts, sequences,
)

# Higher than any mapping quality a read can carry, so every read is rejected.
REJECTS_EVERYTHING = 61


# Both kinds of row have the same shape, so this lists the ones indexed by read
# length. Every other row is indexed by reference position and runs only as far
# as its own reference; a length is not a position, so these rows hold data to
# their full width.
PER_LENGTH = ("reads/lengths",)


def rectangular(output):
    """The arrays with a row per reference, indexed by position or by length."""
    return {k: d[:] for k, d in datasets_of(output).items() if d.ndim == 2}


def counts(output):
    """The rows holding counts, whether a reference's own length bounds them or
    a read length does. These are the ones a zero is meaningful in."""
    return {k: v for k, v in rectangular(output).items() if k not in RATES}


def per_base(output):
    """The arrays a reference's own length bounds, which are the padded ones."""
    return {k: v for k, v in rectangular(output).items() if k not in PER_LENGTH}


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
    return width if field in PER_LENGTH else ref_len


def per_reference(output):
    """The arrays with one value for each reference."""
    return {k: d[:] for k, d in datasets_of(output).items() if d.ndim == 1}


@pytest.fixture
def ragged(datasets, tmp_path):
    """Lengths ranging sixtyfold, so most rows are mostly padding."""
    data = datasets("ragged")
    output = tmp_path / "ragged.h5"
    run_cmuts(data, output)
    return data, output


# ---------------------------------------------------------------------------
# A reference shorter than the row it is written to
# ---------------------------------------------------------------------------


def test_positions_past_a_reference_are_nan(ragged):
    data, output = ragged
    lengths = {name: len(seq) for name, seq in sequences(data.fasta).items()}
    row_of = rows_by_name(data.fasta)

    with h5py.File(output, "r") as handle:
        for field, values in per_base(handle).items():
            width = values.shape[1]
            shorter = [n for n, ln in lengths.items() if ln < width]
            assert shorter, f"{field}: no reference is narrower than the widest"

            for name in shorter:
                tail = values[row_of[name]][lengths[name]:]
                assert np.isnan(tail).all(), f"{field}: {name} is padded with {tail[:4]}"


def test_positions_within_a_reference_are_never_nan(ragged):
    data, output = ragged
    lengths = {name: len(seq) for name, seq in sequences(data.fasta).items()}
    reached = references_with_reads(data.bam)
    row_of = rows_by_name(data.fasta)
    assert reached, "the dataset under test has no reference any read reached"

    with h5py.File(output, "r") as handle:
        for field, values in padded(handle).items():
            for name in reached:
                within = values[row_of[name]][:lengths[name]]
                assert not np.isnan(within).any(), f"{field}: {name} has a hole in it"


# ---------------------------------------------------------------------------
# A reference no read was aligned to
# ---------------------------------------------------------------------------


def test_a_reference_with_no_reads_is_zero_over_its_own_bases(datasets, tmp_path):
    data = datasets("patchy")
    output = tmp_path / "patchy.h5"
    run_cmuts(data, output)

    lengths = {name: len(seq) for name, seq in sequences(data.fasta).items()}
    widest = max(lengths.values())
    reached = references_with_reads(data.bam)
    row_of = rows_by_name(data.fasta)
    missing = [name for name in row_of if name not in reached]

    assert missing, "the dataset under test covers every reference"
    assert any(lengths[name] < widest for name in missing), \
        "no uncovered reference is short enough to carry padding"

    with h5py.File(output, "r") as handle:
        for field, values in counts(handle).items():
            width = values.shape[1]
            for name in missing:
                extent = row_extent(field, lengths[name], width)

                assert (values[row_of[name]][:extent] == 0).all(), \
                    f"{field}: {name} has no reads and is not zero"

        for field, values in padded(handle).items():
            for name in missing:
                assert np.isnan(values[row_of[name]][lengths[name]:]).all(), \
                    f"{field}: {name} has padding that is not NaN"

        # A count of reads is zero where no read arrived: NaN would say nothing
        # that zero does not.
        for field, values in per_reference(handle).items():
            for name in missing:
                assert values[row_of[name]] == 0, \
                    f"{field}: {name} has no reads but is not zero"


def test_an_uncovered_reference_of_full_length_holds_no_nan(datasets, tmp_path):
    data = datasets("flat")
    output = tmp_path / "flat.h5"
    run_cmuts(data, output)

    reached = references_with_reads(data.bam)
    row_of = rows_by_name(data.fasta)
    missing = [name for name in row_of if name not in reached]

    assert missing, "the dataset under test covers every reference"

    with h5py.File(output, "r") as handle:
        for field, values in padded(handle).items():
            for name in missing:
                assert not np.isnan(values[row_of[name]]).any(), \
                    f"{field}: {name} holds a NaN"

        for field, values in counts(handle).items():
            for name in missing:
                assert (values[row_of[name]] == 0).all(), f"{field}: {name} is not zero"


def test_a_reference_whose_reads_were_all_rejected_is_zero(datasets, tmp_path):
    data = datasets("ragged")
    output = tmp_path / "rejected.h5"
    summary = run_cmuts(data, output, min_mapq=REJECTS_EVERYTHING)

    assert summary.kept == 0, "the filter under test let something through"

    lengths = {name: len(seq) for name, seq in sequences(data.fasta).items()}
    reached = references_with_reads(data.bam)
    row_of = rows_by_name(data.fasta)
    assert reached, "the dataset under test has no reference any read reached"

    with h5py.File(output, "r") as handle:
        for field, values in padded(handle).items():
            for name in reached:
                within = values[row_of[name]][:lengths[name]]
                assert not np.isnan(within).any(), f"{field}: {name} went to NaN"

        for field, values in counts(handle).items():
            width = values.shape[1]
            for name in reached:
                within = values[row_of[name]][:row_extent(field, lengths[name], width)]
                assert (within == 0).all(), f"{field}: {name} counted something"

        for field, values in per_reference(handle).items():
            for name in reached:
                assert values[row_of[name]] == 0 or field == "reads/rejected", \
                    f"{field}: {name} is not zero"


# ---------------------------------------------------------------------------
# What NaN means in a rate
# ---------------------------------------------------------------------------


def test_raising_the_minimum_depth_only_adds_nan_to_a_rate(datasets, tmp_path):
    """Separates the two reasons a rate is NaN by the coverage, which is NaN
    outside a reference only."""
    data = datasets("patchy")
    seen = {}

    for depth in (0, 1, 5):
        output = tmp_path / f"depth{depth}.h5"
        run_cmuts(data, output, min_depth=depth)

        with h5py.File(output, "r") as handle:
            reactivity = handle["reactivity"][:]
            error = handle["error"][:]
            coverage = handle["coverage"][:]

        assert np.array_equal(np.isnan(reactivity), np.isnan(error)), \
            f"depth {depth}: the rate and its error disagree about what is known"

        # Padding is not part of any reference, so no rate is had there either.
        assert np.isnan(reactivity[np.isnan(coverage)]).all(), \
            f"depth {depth}: a rate outside a reference is not NaN"

        finite = reactivity[~np.isnan(reactivity)]
        assert finite.size, f"depth {depth}: no position carries a rate at all"
        assert ((finite >= 0) & (finite <= 1)).all(), \
            f"depth {depth}: a rate lies outside zero to one"

        seen[depth] = np.isnan(reactivity)

    assert (seen[1] >= seen[0]).all(), "a higher depth recovered a rate"
    assert (seen[5] >= seen[1]).all(), "a higher depth recovered a rate"
    assert seen[5].sum() > seen[0].sum(), "the depths under test discard the same rates"


def test_an_error_at_a_full_read_of_depth_is_at_most_a_half(datasets, tmp_path):
    """Runs at --min-depth 1, below which the standard error of a proportion
    can exceed a half."""
    data = datasets("patchy")
    output = tmp_path / "bounded.h5"
    run_cmuts(data, output, min_depth=1)

    with h5py.File(output, "r") as handle:
        error = handle["error"][:]

    finite = error[~np.isnan(error)]

    assert finite.size, "no position carries an error at all"
    assert finite.max() <= 0.5, f"an error of {finite.max()} over a whole read"
