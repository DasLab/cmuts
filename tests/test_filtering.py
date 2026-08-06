"""Filtering agrees with samtools, and every read is accounted for."""

import pytest

from conftest import SHAPES
from support import run_cmuts, samtools_kept

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
