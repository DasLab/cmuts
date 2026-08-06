"""Filtering agrees with samtools, and every read is accounted for."""

import contextlib
import os

import pytest

from conftest import SHAPES
from support import Dataset, converted, outputs_agree, run_cmuts, samtools_kept


@contextlib.contextmanager
def unreachable(path):
    """Hides a file for the duration, so nothing can quietly fall back to it."""
    aside = path.with_suffix(path.suffix + ".aside")
    os.rename(path, aside)
    try:
        yield
    finally:
        os.rename(aside, path)

# htslib reads all three, and which one a file is in cannot change what
# survives a filter.
FORMATS = ["bam", "sam", "cram"]

# Combinations chosen so that each criterion is exercised alone, at a boundary,
# and alongside the others.
FILTERS = [
    {},
    {"min_mapq": 1},
    {"min_mapq": 30},
    {"min_mapq": 60},
    {"strand": "forward"},
    {"strand": "reverse"},
    {"min_length": 200},
    {"max_length": 300},
    {"min_length": 150, "max_length": 400},
    {"min_length": 9000},                      # admits nothing
    {"min_mapq": 30, "strand": "reverse", "min_length": 100, "max_length": 500},
]


def describe(filters):
    return ",".join(f"{k}={v}" for k, v in filters.items()) or "unfiltered"


@pytest.mark.parametrize("filters", FILTERS, ids=describe)
@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_matches_samtools(datasets, tmp_path, shape, filters):
    data = datasets(shape)
    summary = run_cmuts(data, tmp_path / "out.h5", **filters)

    assert summary.kept == samtools_kept(data, **filters), "surviving reads"

    # Every mapped read is either kept or rejected; nothing goes missing and
    # nothing is counted twice.
    assert summary.kept + summary.rejected == data.mapped, "reads accounted for"

    # Unmapped reads align nowhere, so no filter can touch them.
    assert summary.unmapped == data.unmapped, "unmapped reads"

    # A reference that received any mapped read gets a row, whether or not
    # anything survived the filter.
    assert summary.rows == data.touched, "references written"


@pytest.mark.parametrize("shape", ["plain", "lowqual"])
def test_rejecting_everything_leaves_a_valid_file(datasets, tmp_path, shape):
    data = datasets(shape)
    summary = run_cmuts(data, tmp_path / "out.h5", min_length=9000)

    assert summary.kept == 0
    assert summary.rejected == data.mapped
    assert summary.rows == data.touched


@pytest.mark.parametrize("filters", FILTERS, ids=describe)
@pytest.mark.parametrize("fmt", FORMATS)
def test_every_format_gives_the_same_answer(datasets, tmp_path, fmt, filters):
    """Filtering is over what a record says, not how it was stored."""
    plain = datasets("plain")
    data = converted(plain, tmp_path, fmt)
    summary = run_cmuts(data, tmp_path / "out.h5", **filters)

    assert summary.kept == samtools_kept(data, **filters)
    assert summary.kept + summary.rejected == data.mapped
    assert summary.unmapped == data.unmapped
    assert summary.rows == data.touched

    run_cmuts(plain, tmp_path / "bam.h5", **filters)
    assert outputs_agree(tmp_path / "out.h5", tmp_path / "bam.h5")


def test_cram_decodes_against_the_reference_it_was_given(datasets, tmp_path):
    """A CRAM names its reference in its own header, by a path recorded when it
    was written and a checksum that may be looked up remotely. Neither need be
    the reference asked for. Here the recorded path is made unreachable, so
    only --fasta can answer."""
    data = converted(datasets("plain"), tmp_path, "cram")

    moved = tmp_path / "elsewhere.fasta"
    moved.write_bytes(data.fasta.read_bytes())
    hidden = Dataset(bam=data.bam, fasta=moved, mapped=data.mapped,
                     unmapped=data.unmapped, touched=data.touched)

    expected = run_cmuts(datasets("plain"), tmp_path / "bam.h5")

    with unreachable(data.fasta):
        assert run_cmuts(hidden, tmp_path / "cram.h5") == expected
