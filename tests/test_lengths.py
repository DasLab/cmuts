"""The read-length histogram, against what samtools makes of the same file.

There is one bin per stored length, from one to twice the longest reference,
and the same bins in every row. A read longer than the last bin is counted in
no bin, leaving its row summing to fewer than the reads counted for that
reference. Only the reads that passed the filter are counted, so samtools is
given the same criteria.
"""

from collections import namedtuple

import numpy as np
import pytest

from oracle import reference_lengths, rows_by_name, samtools_length_histogram
from outputs import COUNTED, LENGTHS, field_of
from programs import run_cmuts

# Unfiltered, and above the mapping quality of most reads the generator
# writes, so the histogram is checked over two different sets of reads.
MAPPING_QUALITIES = (0, 30)

# What a comparison against samtools covered: the references it ran over, and
# the reads too long for any bin to hold.
Compared = namedtuple("Compared", "references outside")


def expected_row(histogram: dict, width: int) -> np.ndarray:
    """Builds the row samtools implies. Bin i holds length i + 1, and a length
    with no bin is dropped."""
    row = np.zeros(width)

    for length, count in histogram.items():
        if 1 <= length <= width:
            row[length - 1] += count

    return row


def compare_against_samtools(output, data, min_mapq: int) -> Compared:
    """Asserts every reference's row against samtools. The criterion is passed
    to both, since the two do not share one."""
    expected = samtools_length_histogram(data, min_mapq=min_mapq)
    lengths = reference_lengths(data.fasta)
    row_of = rows_by_name(data.fasta)
    outside = 0

    written = field_of(output, LENGTHS)
    width = written.shape[1]

    assert width == 2 * max(lengths.values())

    for name, histogram in expected.items():
        assert np.array_equal(written[row_of[name]], expected_row(histogram, width)), name

        outside += sum(count for length, count in histogram.items() if length > width)

    return Compared(references=len(expected), outside=outside)


@pytest.mark.parametrize("min_mapq", MAPPING_QUALITIES)
def test_histogram_matches_samtools(data, falsifiable, tmp_path, min_mapq):
    output = tmp_path / "out.h5"
    run_cmuts(data, output, min_mapq=min_mapq)

    falsifiable(compare_against_samtools(output, data, min_mapq).references > 0)


def test_a_read_longer_than_the_range_is_counted_only_by_the_total(data, falsifiable,
                                                                   tmp_path):
    """A row is short by exactly the reads that no bin could hold. Where no
    read is too long, every row sums to the reads counted for its reference.
    """
    output = tmp_path / "out.h5"
    run_cmuts(data, output, min_mapq=0)

    outside = compare_against_samtools(output, data, min_mapq=0).outside

    falsifiable(outside > 0)

    # Cast to signed, so that a row holding more reads than were counted gives
    # a negative difference. The unsigned difference would wrap to the top of
    # the range.
    counted = field_of(output, COUNTED).astype(np.int64)
    written = field_of(output, LENGTHS).sum(axis=1).astype(np.int64)

    missing = counted - written

    assert (missing >= 0).all()
    assert missing.sum() == outside
