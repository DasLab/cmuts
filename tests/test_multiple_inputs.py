"""Several files counted as one.

Each reference has a single row in the output no matter how many files its reads
came from, so cmuts-hmm merges the files on the reference instead of running them
one after another, which would overwrite the first row with the second.
"""

from dataclasses import replace

import h5py
import numpy as np
import pytest

from datasets import DATASETS
from support import (
    counted_fields, dealt_out, header_of, outputs_agree, reheadered, run_cmuts,
    try_cmuts,
)

# Two and three leave every file holding a share of every busy reference;
# sixteen leaves most files holding nothing for most references.
PARTS = [2, 3, 8, 16]


def values(path, field):
    with h5py.File(path, "r") as output:
        return output[field][:]


# ---------------------------------------------------------------------------
# Splitting the same alignments
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("parts", PARTS)
@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_split_file_counts_the_same_as_the_whole(datasets, falsifiable, tmp_path,
                                                   name, parts):
    data = datasets(name)
    whole = run_cmuts(data, tmp_path / "whole.h5")
    split = run_cmuts(dealt_out(data, tmp_path, parts), tmp_path / "split.h5")

    falsifiable(whole.kept > 0)

    assert split == whole
    assert outputs_agree(tmp_path / "whole.h5", tmp_path / "split.h5")


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_the_order_of_the_input_files_does_not_matter(datasets, falsifiable, tmp_path,
                                                      name):
    data = datasets(name)
    split = dealt_out(data, tmp_path, 4)

    forwards = run_cmuts(split, tmp_path / "forwards.h5")
    run_cmuts(replace(split, bams=tuple(reversed(split.bams))), tmp_path / "backwards.h5")

    falsifiable(forwards.kept > 0)

    assert outputs_agree(tmp_path / "forwards.h5", tmp_path / "backwards.h5")


# ---------------------------------------------------------------------------
# Repeating them
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_the_same_file_twice_counts_everything_twice(datasets, falsifiable, tmp_path,
                                                     name):
    data = datasets(name)
    once = run_cmuts(data, tmp_path / "once.h5")
    twice = run_cmuts(replace(data, bams=data.bams * 2), tmp_path / "twice.h5")

    # Where the first run counted nothing, both totals are zero whether or not
    # cmuts-hmm read the second file.
    falsifiable(once.kept > 0)

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


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_file_declaring_different_references_is_refused(datasets, falsifiable,
                                                          tmp_path, name):
    data = datasets(name)
    renamed = reheadered(data, tmp_path,
                         lambda text: text.replace("SN:ref0000", "SN:other000"))

    # The two files declare different references only where the dataset has a
    # reference the rename matched.
    falsifiable(header_of(renamed) != header_of(data))

    attempt = try_cmuts(replace(data, bams=data.bams + renamed.bams),
                        tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "must be in the same order" in attempt.stderr


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_an_unsorted_file_is_reported_by_path(datasets, falsifiable, tmp_path, name):
    data = datasets(name)
    unsorted = reheadered(data, tmp_path,
                          lambda text: text.replace("SO:coordinate", "SO:unknown"))

    falsifiable(header_of(unsorted) != header_of(data))

    attempt = try_cmuts(replace(data, bams=data.bams + unsorted.bams),
                        tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert str(unsorted.bam) in attempt.stderr
    assert "not coordinate sorted" in attempt.stderr
