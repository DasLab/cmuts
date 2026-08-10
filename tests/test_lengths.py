"""The read-length histogram, against what samtools makes of the same file.

One bin per stored length from zero to twice the reference, and a last bin for
anything longer. What is counted is what survived the filter, so the oracle is
given the criteria cmuts was given and measures the sequence column itself.
"""

import h5py
import numpy as np
import pytest

from conftest import SHAPES
from support import (
    rows_by_name, run_cmuts, samtools_length_histogram, sequences,
)


def expected_row(histogram, ref_len, width):
    """The row samtools implies: every length in its own bin, and everything
    past the range gathered into the last one."""
    row = np.zeros(width)
    overflow = 2 * ref_len + 1

    for length, count in histogram.items():
        row[min(length, overflow)] += count

    return row


def compare(output, data, min_mapq):
    """Every reference's row against the histogram samtools implies, under the
    same criterion cmuts was given. The criterion is passed to both rather than
    left to a default, the two not sharing one."""
    expected = samtools_length_histogram(data, min_mapq=min_mapq)
    lengths = {name: len(seq) for name, seq in sequences(data.fasta).items()}
    overflowed = 0

    with h5py.File(output, "r") as handle:
        row_of = rows_by_name(handle)
        written = handle["read_lengths"][:]

        for name, histogram in expected.items():
            ref_len = lengths[name]
            extent = 2 * ref_len + 2
            row = written[row_of[name]]

            assert np.array_equal(row[:extent], expected_row(histogram, ref_len, extent)), \
                f"{name}: histogram disagrees with samtools"

            overflowed += sum(n for n in histogram if n > 2 * ref_len + 1)

    return overflowed


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


def test_the_histogram_sums_to_the_reads_counted(datasets, tmp_path):
    data = datasets("ragged")
    output = tmp_path / "sums.h5"
    run_cmuts(data, output, min_mapq=0)

    with h5py.File(output, "r") as handle:
        reads = handle["reads"][:]
        written = handle["read_lengths"][:]
        reached = ~np.isnan(reads)

        assert np.array_equal(np.nansum(written[reached], axis=1), reads[reached])


def test_a_read_longer_than_the_range_lands_in_the_overflow_bin(datasets, tmp_path):
    """Soft-clipped bases are stored but align nowhere, so a read can be far
    longer than twice the reference it was placed on. The last bin is where
    those go, and nothing past it exists to hold them."""
    data = datasets("overflowing")
    output = tmp_path / "overflowing.h5"
    run_cmuts(data, output, min_mapq=0)

    overflowed = compare(output, data, min_mapq=0)

    assert overflowed, "the shape under test produced no read past the range"
