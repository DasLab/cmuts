"""Filtering agrees with samtools, and every read is accounted for."""

import pytest

from datasets import DATASETS
from filters import FILTERS, UNFILTERED, criteria, describe_filters
from support import (
    MAPQ_COLUMN,
    UNAVAILABLE_MAPQ,
    assert_counts_agree,
    records,
    run_cmuts,
    samtools_kept,
    with_secondary,
)


@pytest.mark.parametrize("filters", FILTERS, ids=describe_filters)
@pytest.mark.parametrize("name", sorted(DATASETS))
def test_read_counts_match_samtools(datasets, falsifiable, tmp_path, name, filters):
    data = datasets(name)

    falsifiable(data.mapped > 0)

    summary = run_cmuts(data, tmp_path / "out.h5", **criteria(filters))

    assert_counts_agree(summary, data, criteria(filters))


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_an_unavailable_mapping_quality_is_refused(datasets, falsifiable, tmp_path, name):
    """Counted from the records here rather than taken from samtools_kept,
    which is told of the same divergence and would agree with a filter that
    never applied it."""
    data = datasets(name)
    mapped = records(data.bam, "-F", "0x104")
    unavailable = [line for line in mapped
                   if int(line.split("\t")[MAPQ_COLUMN]) == UNAVAILABLE_MAPQ]

    falsifiable(len(unavailable) > 0)

    summary = run_cmuts(data, tmp_path / "out.h5", **UNFILTERED)

    assert summary.kept == len(mapped) - len(unavailable), \
        "a read of unavailable mapping quality survived a threshold of zero"


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_rejecting_everything_leaves_a_valid_file(datasets, falsifiable, tmp_path, name):
    data = datasets(name)

    falsifiable(data.mapped > 0)

    summary = run_cmuts(data, tmp_path / "out.h5", min_length=9000)

    assert summary.kept == 0
    assert summary.rejected == data.mapped
    assert summary.rows == data.touched


# One in three of the mapped reads, so that a dataset holding any leaves both a
# marked read and an unmarked one.
EVERY_THIRD = 3


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_secondary_alignments_are_refused(datasets, falsifiable, tmp_path, name):
    data, marked = with_secondary(datasets(name), tmp_path, every=EVERY_THIRD)

    falsifiable(marked > 0)

    summary = run_cmuts(data, tmp_path / "out.h5", **UNFILTERED)

    assert summary.kept == samtools_kept(data, **UNFILTERED), "agreement with samtools"

    # They are refused rather than overlooked, so they show up as rejected.
    # Others may be refused alongside them, which is why the count is a bound
    # and not the total: a read of unavailable mapping quality is refused for
    # a reason of its own.
    assert summary.rejected >= marked, "a marked read was not refused"
    assert summary.kept + summary.rejected == data.mapped, "reads accounted for"
