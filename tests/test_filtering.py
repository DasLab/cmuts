"""Filtering agrees with samtools, and every read is accounted for."""

import pytest

from alignments import mark_secondary
from filters import FILTERS, UNFILTERED, criteria, describe_filters
from oracle import (
    NOT_COUNTED_FLAGS,
    UNAVAILABLE_MAPQ,
    assert_counts_agree,
    records,
    samtools_kept,
)
from outputs import read_summary
from programs import run_cmuts


@pytest.mark.parametrize("filters", FILTERS, ids=describe_filters)
def test_read_counts_match_samtools(data, falsifiable, tmp_path, filters):
    falsifiable(data.mapped > 0)

    summary = read_summary(run_cmuts(data, tmp_path / "out.h5", **criteria(filters)))

    assert_counts_agree(summary, data, criteria(filters))


def test_an_unavailable_mapping_quality_is_refused(data, falsifiable, tmp_path):
    """Counts the affected reads from the records directly. samtools_kept
    applies the same divergence, so it would agree with a filter that never
    applied it."""
    counted = records(data.bam, *NOT_COUNTED_FLAGS)
    unavailable = [record for record in counted if record.mapq == UNAVAILABLE_MAPQ]

    falsifiable(len(unavailable) > 0)

    summary = read_summary(run_cmuts(data, tmp_path / "out.h5", **UNFILTERED))

    assert summary.kept == len(counted) - len(unavailable)


def test_rejecting_everything_leaves_a_valid_file(data, falsifiable, tmp_path):
    falsifiable(data.mapped > 0)

    summary = read_summary(run_cmuts(data, tmp_path / "out.h5", min_length=9000))

    assert summary.kept == 0
    assert summary.rejected == data.mapped
    assert summary.rows == data.touched


# One in three of the mapped reads, so that a dataset holding any leaves both a
# marked read and an unmarked one.
EVERY_THIRD = 3


def test_secondary_alignments_are_refused(data, falsifiable, tmp_path):
    marked_data, marked = mark_secondary(data, tmp_path, every=EVERY_THIRD)

    falsifiable(marked > 0)

    summary = read_summary(run_cmuts(marked_data, tmp_path / "out.h5", **UNFILTERED))

    assert summary.kept == samtools_kept(marked_data, **UNFILTERED)

    # A marked read is rejected and not skipped, so it reaches the rejected
    # total. Other reads may be rejected alongside it, so the count of marked
    # reads is a lower bound and not the total.
    assert summary.rejected >= marked
    assert summary.kept + summary.rejected == marked_data.mapped
