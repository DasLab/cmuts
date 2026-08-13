"""Filtering agrees with samtools, and every read is accounted for."""

import pytest

from datasets import DATASETS
from filters import FILTERS, UNFILTERED, criteria, describe_filters
from support import (
    assert_counts_agree,
    generate,
    run_cmuts,
    samtools_kept,
    with_secondary,
)


@pytest.mark.parametrize("filters", FILTERS, ids=describe_filters)
@pytest.mark.parametrize("name", sorted(DATASETS))
def test_read_counts_match_samtools(datasets, tmp_path, name, filters):
    data = datasets(name)
    summary = run_cmuts(data, tmp_path / "out.h5", **criteria(filters))

    assert_counts_agree(summary, data, criteria(filters))


def test_an_unavailable_mapping_quality_is_refused(tmp_path):
    """Asserted against the summary and not against samtools, which admits
    MAPQ 255 at every threshold."""
    data = generate(tmp_path, "unavailable", seed=112, references=8,
                    reads_per_ref=10, mapq=255, unmapped=0)

    summary = run_cmuts(data, tmp_path / "out.h5", **UNFILTERED)

    assert data.mapped > 0
    assert summary.kept == 0
    assert summary.rejected == data.mapped


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_rejecting_everything_leaves_a_valid_file(datasets, tmp_path, name):
    data = datasets(name)
    summary = run_cmuts(data, tmp_path / "out.h5", min_length=9000)

    assert summary.kept == 0
    assert summary.rejected == data.mapped
    assert summary.rows == data.touched


# One in three of the mapped reads, so that a dataset holding any leaves both a
# marked read and an unmarked one.
EVERY_THIRD = 3


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_secondary_alignments_are_refused(datasets, checked, tmp_path, name):
    data, marked = with_secondary(datasets(name), tmp_path, every=EVERY_THIRD)

    checked(marked)

    summary = run_cmuts(data, tmp_path / "out.h5", **UNFILTERED)

    assert summary.kept == data.mapped - marked, "surviving reads"
    assert summary.kept == samtools_kept(data, **UNFILTERED), "agreement with samtools"

    # They are refused rather than overlooked, so they show up as rejected.
    assert summary.kept + summary.rejected == data.mapped, "reads accounted for"
