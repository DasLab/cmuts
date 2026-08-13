"""The read-length histogram, against what samtools makes of the same file.

One bin per stored length from one to twice the longest reference, the same
bins in every row: a length is not a position, so a column means one length in
every row alike. A read longer than the last bin is counted in none, leaving a
row summing to fewer than the reads total. Only what survived the filter is
counted, so samtools is given the same criteria.
"""

from collections import namedtuple

import h5py
import numpy as np
import pytest

from datasets import DATASETS
from support import (
    reference_lengths, rows_by_name, run_cmuts, samtools_length_histogram,
)

# Unfiltered, and above the mapping quality most of what the generator writes
# carries, so the histogram is checked over two different sets of reads.
CRITERIA = (0, 30)

# What a comparison against samtools covered: the references it ran over, and
# the reads too long for any bin to hold.
Compared = namedtuple("Compared", "references outside")


def expected_row(histogram, width):
    """The row samtools implies. Bin i holds length i + 1, the histogram
    beginning at 1 rather than 0, and a length with no bin is dropped."""
    row = np.zeros(width)

    for length, count in histogram.items():
        if 1 <= length <= width:
            row[length - 1] += count

    return row


def compare(output, data, min_mapq) -> Compared:
    """Checks every reference's row against samtools. The criterion is passed
    to both rather than left to a default, the two not sharing one."""
    expected = samtools_length_histogram(data, min_mapq=min_mapq)
    lengths = reference_lengths(data.fasta)
    row_of = rows_by_name(data.fasta)
    outside = 0

    with h5py.File(output, "r") as handle:
        written = handle["reads/lengths"][:]
        width = written.shape[1]

        assert width == 2 * max(lengths.values()), "the row is not the widest reference"

        for name, histogram in expected.items():
            row = written[row_of[name]]

            assert np.array_equal(row, expected_row(histogram, width)), \
                f"{name}: histogram disagrees with samtools"

            outside += sum(c for n, c in histogram.items() if n > width)

    return Compared(references=len(expected), outside=outside)


@pytest.mark.parametrize("min_mapq", CRITERIA)
@pytest.mark.parametrize("name", sorted(DATASETS))
def test_histogram_matches_samtools(datasets, falsifiable, tmp_path, name, min_mapq):
    data = datasets(name)
    output = tmp_path / f"{name}.h5"
    run_cmuts(data, output, min_mapq=min_mapq)

    falsifiable(compare(output, data, min_mapq).references > 0)


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_read_longer_than_the_range_is_counted_only_by_the_total(datasets, falsifiable,
                                                                   tmp_path, name):
    """A row is short by the reads no bin could hold and by nothing else, which
    where none is too long leaves every row summing to the reads counted for
    its reference.
    """
    data = datasets(name)
    output = tmp_path / f"{name}.h5"
    run_cmuts(data, output, min_mapq=0)

    outside = compare(output, data, min_mapq=0).outside

    falsifiable(outside > 0)

    with h5py.File(output, "r") as handle:
        # Signed, so that a row holding more reads than were counted comes out
        # negative instead of wrapping to the top of the unsigned range.
        counted = handle["reads/counted"][:].astype(np.int64)
        written = handle["reads/lengths"][:].sum(axis=1).astype(np.int64)

    missing = counted - written

    assert (missing >= 0).all(), "a row holds more reads than were counted"
    assert missing.sum() == outside, \
        "the reads outside the range are not what the total is short by"
