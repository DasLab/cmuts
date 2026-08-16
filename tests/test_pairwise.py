"""What a pairwise count writes, whatever the alignment it was given.

`--pairwise` adds a value for every ordered pair of a reference's positions: how far the
two are modified together, and the reads behind that. Each test here asserts something
that holds of any alignment, so none of them looks at what the correlations came to.

The output goes as the square of the longest reference, so these run over the datasets
with short ones. What is counted does not depend on the format the alignment arrives in,
so one format is enough.
"""

import h5py
import numpy as np
import pytest

from alignments import NATIVE
from outputs import ALL_FIELDS, COUNTED, arrays_agree, field_of

from programs import run_cmuts

CORRELATION = "pairwise/correlation"
COVERAGE = "pairwise/coverage"

# Datasets whose references are short enough to square.
SHORT = ("tiny", "deep")

# More reads than any dataset here holds, so no pair meets it.
BEYOND_EVERY_DEPTH = 10_000_000


@pytest.fixture(params=SHORT)
def short(request, catalogue):
    """Yields a dataset short enough for a square of it to stay small."""
    return catalogue(request.param, NATIVE)


@pytest.fixture
def output(short, tmp_path):
    """A pairwise count of that dataset, as the path it was written to."""
    return run_cmuts(short, tmp_path / "pairwise.h5", pairwise=True)


def holds(output, name) -> bool:
    """Returns whether an output carries a dataset."""
    with h5py.File(output, "r") as handle:
        return name in handle


def mirrored(square) -> np.ndarray:
    """Returns the square with each reference's pairs read the other way round."""
    return square.transpose(0, 2, 1)


# ---------------------------------------------------------------------------
# What is written
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", [CORRELATION, COVERAGE])
def test_a_pairwise_field_is_absent_unless_asked_for(short, falsifiable, tmp_path, name):
    plain = run_cmuts(short, tmp_path / "plain.h5")

    falsifiable(True)

    assert not holds(plain, name)


def test_a_square_is_written_for_every_reference(output, falsifiable):
    """As wide as the rows every other field writes, so a column means the same
    position throughout the output."""
    columns = field_of(output, "coverage").shape[1]
    correlation = field_of(output, CORRELATION)

    falsifiable(columns > 0)

    assert correlation.shape == (correlation.shape[0], columns, columns)


@pytest.mark.parametrize("name", ALL_FIELDS)
def test_asking_for_pairs_leaves_the_other_fields_alone(short, falsifiable, tmp_path,
                                                        name):
    plain = run_cmuts(short, tmp_path / "plain.h5")
    paired = run_cmuts(short, tmp_path / "paired.h5", pairwise=True)

    falsifiable(True)

    assert arrays_agree(field_of(plain, name), field_of(paired, name))


# ---------------------------------------------------------------------------
# What holds of the squares
# ---------------------------------------------------------------------------


def test_a_correlation_reads_the_same_either_way_round(output, falsifiable):
    correlation = field_of(output, CORRELATION)

    falsifiable(correlation.size > 0)

    assert np.array_equal(correlation, mirrored(correlation), equal_nan=True)


def test_a_pair_is_covered_the_same_either_way_round(output, falsifiable):
    coverage = field_of(output, COVERAGE)

    falsifiable(coverage.size > 0)

    assert np.array_equal(coverage, mirrored(coverage))


def test_no_correlation_lies_outside_the_range(output, falsifiable):
    correlation = field_of(output, CORRELATION)
    known = np.isfinite(correlation)

    falsifiable(known.any())

    assert np.all(np.abs(correlation[known]) <= 1.0)


def test_a_pair_is_covered_by_no_more_reads_than_were_counted(output, falsifiable):
    """The reads reaching two positions are among those the run counted at all."""
    coverage = field_of(output, COVERAGE)
    counted = field_of(output, COUNTED)

    falsifiable(counted.sum() > 0)

    assert np.all(coverage.max(axis=(1, 2)) <= counted)


def test_a_pair_below_the_depth_carries_no_correlation(short, falsifiable, tmp_path):
    output = run_cmuts(short, tmp_path / "deep.h5", pairwise=True,
                       min_depth=BEYOND_EVERY_DEPTH)

    falsifiable(True)

    assert np.all(np.isnan(field_of(output, CORRELATION)))
