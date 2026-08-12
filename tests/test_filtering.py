"""Filtering agrees with samtools, and every read is accounted for."""

import contextlib
import os
from dataclasses import replace

import pytest

from conftest import SHAPES
from support import (
    assert_counts_agree,
    converted,
    outputs_agree,
    run_cmuts,
    samtools_kept,
    with_secondary,
)


@contextlib.contextmanager
def unreachable(path):
    """Renames a file aside for the duration, so nothing can fall back to it."""
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


# Criteria are given in full at every call, so no default is relied on here.
UNFILTERED = {"min_mapq": 0, "strand": "both", "min_length": 0, "max_length": 0}


def criteria(filters):
    return {**UNFILTERED, **filters}


def describe(filters):
    return ",".join(f"{k}={v}" for k, v in filters.items()) or "unfiltered"


@pytest.mark.parametrize("filters", FILTERS, ids=describe)
@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_read_counts_match_samtools(datasets, tmp_path, shape, filters):
    data = datasets(shape)
    summary = run_cmuts(data, tmp_path / "out.h5", **criteria(filters))

    assert_counts_agree(summary, data, criteria(filters))


@pytest.mark.parametrize("shape", ["plain", "lowqual"])
def test_rejecting_everything_leaves_a_valid_file(datasets, tmp_path, shape):
    data = datasets(shape)
    summary = run_cmuts(data, tmp_path / "out.h5", min_length=9000)

    assert summary.kept == 0
    assert summary.rejected == data.mapped
    assert summary.rows == data.touched


def test_secondary_alignments_are_refused(datasets, tmp_path):
    data, marked = with_secondary(datasets("plain"), tmp_path, every=3)
    assert marked > 0, "the case being tested has to appear in the data"

    summary = run_cmuts(data, tmp_path / "out.h5", **UNFILTERED)

    assert summary.kept == data.mapped - marked, "surviving reads"
    assert summary.kept == samtools_kept(data, **UNFILTERED), "agreement with samtools"

    # They are refused rather than overlooked, so they show up as rejected.
    assert summary.kept + summary.rejected == data.mapped, "reads accounted for"


@pytest.mark.parametrize("filters", FILTERS, ids=describe)
@pytest.mark.parametrize("fmt", FORMATS)
def test_every_format_gives_the_same_answer(datasets, tmp_path, fmt, filters):
    plain = datasets("plain")
    data = converted(plain, tmp_path, fmt)
    summary = run_cmuts(data, tmp_path / "out.h5", **criteria(filters))

    assert_counts_agree(summary, data, criteria(filters))

    run_cmuts(plain, tmp_path / "bam.h5", **criteria(filters))
    assert outputs_agree(tmp_path / "out.h5", tmp_path / "bam.h5")


def test_cram_decodes_against_the_reference_it_was_given(datasets, tmp_path):
    """Renames the reference recorded in the CRAM header aside, leaving cmuts
    no source for the bases but --fasta."""
    data = converted(datasets("plain"), tmp_path, "cram")

    moved = tmp_path / "elsewhere.fasta"
    moved.write_bytes(data.fasta.read_bytes())
    hidden = replace(data, fasta=moved)

    run_cmuts(datasets("plain"), tmp_path / "bam.h5")

    with unreachable(data.fasta):
        run_cmuts(hidden, tmp_path / "cram.h5")

    # Read counts agree whichever reference cmuts decoded against, so only the
    # per-base fields show which one it used.
    assert outputs_agree(tmp_path / "cram.h5", tmp_path / "bam.h5")
