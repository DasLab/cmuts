"""Several files counted as one.

Each reference has a single row in the output no matter how many files its reads
came from, so cmuts merges the files on the reference instead of running them
one after another, which would overwrite the first row with the second.
"""

from dataclasses import replace

import h5py
import numpy as np
import pytest

from support import (
    counted_fields, dealt_out, outputs_agree, reheadered, run_cmuts, try_cmuts,
)


def values(path, field):
    with h5py.File(path, "r") as output:
        return output[field][:]


# ---------------------------------------------------------------------------
# Splitting the same alignments
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("parts", [2, 3, 8])
def test_a_split_file_counts_the_same_as_the_whole(data, tmp_path, parts):
    whole = run_cmuts(data, tmp_path / "whole.h5")
    split = run_cmuts(dealt_out(data, tmp_path, parts), tmp_path / "split.h5")

    assert split == whole
    assert outputs_agree(tmp_path / "whole.h5", tmp_path / "split.h5")


def test_more_files_than_references_still_counts_the_same(datasets, tmp_path):
    """Sixteen files over a sparse dataset, so most hold nothing for most
    references."""
    data = datasets("sparse")

    run_cmuts(data, tmp_path / "whole.h5")
    run_cmuts(dealt_out(data, tmp_path, 16), tmp_path / "split.h5")

    assert outputs_agree(tmp_path / "whole.h5", tmp_path / "split.h5")


def test_the_order_the_files_are_named_in_does_not_matter(data, tmp_path):
    split = dealt_out(data, tmp_path, 4)

    run_cmuts(split, tmp_path / "forwards.h5")
    run_cmuts(replace(split, bams=tuple(reversed(split.bams))), tmp_path / "backwards.h5")

    assert outputs_agree(tmp_path / "forwards.h5", tmp_path / "backwards.h5")


# ---------------------------------------------------------------------------
# Repeating them
# ---------------------------------------------------------------------------


def test_the_same_file_twice_counts_everything_twice(data, tmp_path):
    once = run_cmuts(data, tmp_path / "once.h5")
    twice = run_cmuts(replace(data, bams=data.bams * 2), tmp_path / "twice.h5")

    assert twice.kept == 2 * once.kept
    assert twice.unmapped == 2 * once.unmapped
    assert twice.rows == once.rows

    for field in counted_fields(tmp_path / "once.h5"):
        assert np.array_equal(values(tmp_path / "twice.h5", field),
                              2 * values(tmp_path / "once.h5", field),
                              equal_nan=True), field


# ---------------------------------------------------------------------------
# Files that do not belong together
# ---------------------------------------------------------------------------


def test_a_file_naming_its_references_differently_is_refused(data, tmp_path):
    """A record names its reference by an index into its own header, so
    disagreeing headers put reads in the wrong rows."""
    renamed = reheadered(data, tmp_path,
                         lambda text: text.replace("SN:ref0000", "SN:other000"))

    attempt = try_cmuts(replace(data, bams=data.bams + renamed.bams),
                        tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "must be in the same order" in attempt.stderr


def test_a_file_that_is_not_coordinate_sorted_is_refused_by_name(data, tmp_path):
    unsorted =reheadered(data, tmp_path,
                          lambda text: text.replace("SO:coordinate", "SO:unknown"))

    attempt = try_cmuts(replace(data, bams=data.bams + unsorted.bams),
                        tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert str(unsorted.bam) in attempt.stderr
    assert "not coordinate sorted" in attempt.stderr
