"""The read-length histogram, against what samtools makes of the same file.

One bin per stored length from one to twice the longest reference, the same
bins in every row: a length is not a position, so a column means one length in
every row alike. A read longer than the last bin is counted in none, leaving a
row summing to fewer than the reads total. Only what survived the filter is
counted, so samtools is given the same criteria.
"""

import h5py
import numpy as np
import pytest

from conftest import SHAPES
from support import (
    rows_by_name, run_cmuts, samtools_length_histogram, sequences,
)


def expected_row(histogram, width):
    """The row samtools implies. Bin i holds length i + 1, the histogram
    beginning at 1 rather than 0, and a length with no bin is dropped."""
    row = np.zeros(width)

    for length, count in histogram.items():
        if 1 <= length <= width:
            row[length - 1] += count

    return row


def compare(output, data, min_mapq):
    """Checks every reference's row against samtools, returning how many reads
    fell outside the range of bins. The criterion is passed to both rather than
    left to a default, the two not sharing one."""
    expected = samtools_length_histogram(data, min_mapq=min_mapq)
    lengths = {name: len(seq) for name, seq in sequences(data.fasta).items()}
    row_of = rows_by_name(data.fasta)
    outside = 0

    # A criterion that admits nothing leaves nothing to compare, and every
    # assertion below would hold of an output that counted the wrong thing.
    assert expected, "no read survives the criterion under test"

    with h5py.File(output, "r") as handle:
        written = handle["reads/lengths"][:]
        width = written.shape[1]

        assert width == 2 * max(lengths.values()), "the row is not the widest reference"

        for name, histogram in expected.items():
            row = written[row_of[name]]

            assert np.array_equal(row, expected_row(histogram, width)), \
                f"{name}: histogram disagrees with samtools"

            outside += sum(c for n, c in histogram.items() if n > width)

    return outside


@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_histogram_matches_samtools(datasets, tmp_path, shape):
    data = datasets(shape)
    output = tmp_path / f"{shape}.h5"
    run_cmuts(data, output, min_mapq=0)

    compare(output, data, min_mapq=0)


@pytest.mark.parametrize("shape", ["plain", "clipped"])
def test_histogram_matches_samtools_under_a_filter(datasets, tmp_path, shape):
    data = datasets(shape)
    output = tmp_path / f"{shape}-filtered.h5"
    run_cmuts(data, output, min_mapq=30)

    compare(output, data, min_mapq=30)


def test_each_row_sums_to_the_reads_counted_for_its_reference(datasets, tmp_path):
    data = datasets("ragged")
    output = tmp_path / "sums.h5"
    run_cmuts(data, output, min_mapq=0)

    assert compare(output, data, min_mapq=0) == 0, "this shape overflows the range"

    with h5py.File(output, "r") as handle:
        reads = handle["reads/counted"][:]
        written = handle["reads/lengths"][:]
        reached = ~np.isnan(reads)

        assert np.array_equal(np.nansum(written[reached], axis=1), reads[reached])


def test_a_read_longer_than_the_range_is_counted_only_by_the_total(datasets, tmp_path):
    data = datasets("overflowing")
    output = tmp_path / "overflowing.h5"
    run_cmuts(data, output, min_mapq=0)

    outside = compare(output, data, min_mapq=0)
    assert outside, "the shape under test produced no read past the range"

    with h5py.File(output, "r") as handle:
        reads = handle["reads/counted"][:]
        written = handle["reads/lengths"][:]
        reached = ~np.isnan(reads)

        missing = reads[reached] - np.nansum(written[reached], axis=1)

        assert (missing >= 0).all(), "a row holds more reads than were counted"
        assert missing.sum() == outside, \
            "the reads outside the range are not what the total is short by"
