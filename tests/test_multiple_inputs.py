"""Several files counted as one.

Each reference has a single row in the output no matter how many files its
reads came from, so cmuts-hmm merges the files on the reference. Running them
one after another would overwrite the first row with the second.
"""

from dataclasses import replace

import numpy as np
import pytest

from alignments import replace_header, split_across_files
from oracle import header_text
from outputs import COUNT_FIELDS, field_of, outputs_agree
from programs import run_cmuts, try_cmuts

# Two and three leave every file holding a share of every busy reference.
# Sixteen leaves most files holding nothing for most references.
PARTS = [2, 3, 8, 16]

# Enough files that reversing their order is a different order.
SPLIT_PARTS = 4


# ---------------------------------------------------------------------------
# Splitting one file
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("parts", PARTS)
def test_a_split_file_counts_the_same_as_the_whole(data, falsifiable, tmp_path, parts):
    whole = run_cmuts(data, tmp_path / "whole.h5")
    split = run_cmuts(split_across_files(data, tmp_path, parts), tmp_path / "split.h5")

    falsifiable(whole.kept > 0)

    assert split == whole
    assert outputs_agree(tmp_path / "whole.h5", tmp_path / "split.h5")


def test_the_order_of_the_input_files_does_not_matter(data, falsifiable, tmp_path):
    split = split_across_files(data, tmp_path, SPLIT_PARTS)

    forwards = run_cmuts(split, tmp_path / "forwards.h5")
    run_cmuts(replace(split, bams=tuple(reversed(split.bams))),
              tmp_path / "backwards.h5")

    falsifiable(forwards.kept > 0)

    assert outputs_agree(tmp_path / "forwards.h5", tmp_path / "backwards.h5")


# ---------------------------------------------------------------------------
# Repeating one file
# ---------------------------------------------------------------------------


def test_the_same_file_twice_counts_everything_twice(data, falsifiable, tmp_path):
    once = run_cmuts(data, tmp_path / "once.h5")
    twice = run_cmuts(replace(data, bams=data.bams * 2), tmp_path / "twice.h5")

    # Where the first run kept no reads, both totals are zero whether or not
    # cmuts-hmm read the second file.
    falsifiable(once.kept > 0)

    assert twice.kept == 2 * once.kept
    assert twice.unmapped == 2 * once.unmapped
    assert twice.rows == once.rows

    for field in COUNT_FIELDS:
        assert np.array_equal(field_of(tmp_path / "twice.h5", field),
                              2 * field_of(tmp_path / "once.h5", field),
                              equal_nan=True), field


# ---------------------------------------------------------------------------
# Files that do not belong together
# ---------------------------------------------------------------------------


def test_a_file_declaring_different_references_is_refused(data, falsifiable, tmp_path):
    renamed = replace_header(
        data, tmp_path, lambda text: text.replace("SN:ref0000", "SN:other000"))

    # The two files declare different references only where the dataset has a
    # reference whose name the rename matched.
    falsifiable(header_text(renamed.bam) != header_text(data.bam))

    attempt = try_cmuts(replace(data, bams=data.bams + renamed.bams),
                        tmp_path / "out.h5")

    assert attempt.returncode != 0


def test_an_unsorted_file_among_several_is_refused(data, falsifiable, tmp_path):
    unsorted = replace_header(
        data, tmp_path, lambda text: text.replace("SO:coordinate", "SO:unknown"))

    falsifiable(header_text(unsorted.bam) != header_text(data.bam))

    attempt = try_cmuts(replace(data, bams=data.bams + unsorted.bams),
                        tmp_path / "out.h5")

    assert attempt.returncode != 0
