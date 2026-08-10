"""The read-length histogram, against what samtools makes of the same file.

One bin per stored length from zero to twice the longest reference, the same
bins in every row: a length is not a position, so a column means one length
whatever reference the row belongs to. A read longer than the range is counted
in no bin, so a row sums to the reads it holds and the reads total says how
many fell outside. What is counted is what survived the filter, so the oracle
is given the criteria cmuts was given and measures the sequence column itself.
"""

import h5py
import numpy as np
import pytest

from conftest import SHAPES
from support import (
    rows_by_name, run_cmuts, samtools_length_histogram, sequences,
)


def expected_row(histogram, width):
    """The row samtools implies: every length it reports in its own bin, and
    nothing at all for a length the row has no bin for."""
    row = np.zeros(width)

    for length, count in histogram.items():
        if length < width:
            row[length] += count

    return row


def compare(output, data, min_mapq):
    """Every reference's row against the histogram samtools implies, under the
    same criterion cmuts was given. The criterion is passed to both rather than
    left to a default, the two not sharing one."""
    expected = samtools_length_histogram(data, min_mapq=min_mapq)
    lengths = {name: len(seq) for name, seq in sequences(data.fasta).items()}
    outside = 0

    with h5py.File(output, "r") as handle:
        row_of = rows_by_name(handle)
        written = handle["read_lengths"][:]
        width = written.shape[1]

        assert width == 2 * max(lengths.values()) + 1, "the row is not the widest reference"

        for name, histogram in expected.items():
            row = written[row_of[name]]

            assert np.array_equal(row, expected_row(histogram, width)), \
                f"{name}: histogram disagrees with samtools"

            outside += sum(c for n, c in histogram.items() if n >= width)

    return outside


@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_histogram_matches_samtools(datasets, tmp_path, shape):
    data = datasets(shape)
    output = tmp_path / f"{shape}.h5"
    run_cmuts(data, output, min_mapq=0)

    compare(output, data, min_mapq=0)


@pytest.mark.parametrize("shape", ["plain", "clipped"])
def test_histogram_matches_samtools_under_a_filter(datasets, tmp_path, shape):
    """The histogram counts what the tally saw, so a filter that removes reads
    has to remove them from here too."""
    data = datasets(shape)
    output = tmp_path / f"{shape}-filtered.h5"
    run_cmuts(data, output, min_mapq=30)

    compare(output, data, min_mapq=30)


def test_the_histogram_sums_to_the_reads_it_holds(datasets, tmp_path):
    """Every read of a shape whose lengths are all in range is binned, so the
    rows come to the reads total exactly."""
    data = datasets("ragged")
    output = tmp_path / "sums.h5"
    run_cmuts(data, output, min_mapq=0)

    assert compare(output, data, min_mapq=0) == 0, "this shape overflows the range"

    with h5py.File(output, "r") as handle:
        reads = handle["reads"][:]
        written = handle["read_lengths"][:]
        reached = ~np.isnan(reads)

        assert np.array_equal(np.nansum(written[reached], axis=1), reads[reached])


def test_a_read_longer_than_the_range_is_counted_by_the_total_alone(datasets, tmp_path):
    """A read can be longer than twice the longest reference, soft-clipped bases
    being stored while aligning nowhere. No bin holds it, and what says it was
    there is the reads total standing above the row's own sum."""
    data = datasets("overflowing")
    output = tmp_path / "overflowing.h5"
    run_cmuts(data, output, min_mapq=0)

    outside = compare(output, data, min_mapq=0)
    assert outside, "the shape under test produced no read past the range"

    with h5py.File(output, "r") as handle:
        reads = handle["reads"][:]
        written = handle["read_lengths"][:]
        reached = ~np.isnan(reads)

        missing = reads[reached] - np.nansum(written[reached], axis=1)

        assert (missing >= 0).all(), "a row holds more reads than were counted"
        assert missing.sum() == outside, \
            "the reads outside the range are not what the total is short by"
