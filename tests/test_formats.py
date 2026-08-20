"""The same alignments in each format htslib reads.

The format of a file decides only how the bases and the CIGAR are stored, so a
run over one format must answer as a run over another. Every other
test file runs over each format in turn; what is left here is the comparison
between them, and the reference a CRAM reader is given.
"""

import contextlib
import os
from dataclasses import replace

import pytest

from alignments import CRAM, FORMATS, NATIVE, convert_format
from filters import COMPOUND, criteria, describe_filters
from oracle import assert_counts_agree
from outputs import outputs_agree, read_summary
from programs import run_cmuts

# The criteria a format is tried under. The whole matrix belongs to
# test_filtering.py, which runs it over every dataset. A format can only change
# what a filter reads, so what is needed here is one run reading everything a
# filter can read and one run rejecting every read.
CRITERIA = [{}, COMPOUND]


@pytest.fixture(params=[NATIVE])
def fmt(request):
    """Each test here converts the alignments itself, and compares the result
    against the format they were generated in."""
    return request.param


@contextlib.contextmanager
def moved_aside(path):
    """Renames a file for the duration of the block, so that no run can open
    it."""
    aside = path.with_suffix(path.suffix + ".aside")
    os.rename(path, aside)
    try:
        yield
    finally:
        os.rename(aside, path)


@pytest.mark.parametrize("filters", CRITERIA, ids=describe_filters)
@pytest.mark.parametrize("target", FORMATS)
def test_every_format_gives_the_same_answer(data, falsifiable, tmp_path, target,
                                            filters):
    """The target is named separately from the fmt fixture, which holds the
    alignments in the format they were generated in for the comparison."""
    converted = convert_format(data, tmp_path, target)
    summary = read_summary(
        run_cmuts(converted, tmp_path / "out.h5", **criteria(filters)))

    # Only a run that kept reads compares what the format stored.
    falsifiable(summary.kept > 0)

    assert_counts_agree(summary, converted, criteria(filters))

    run_cmuts(data, tmp_path / "bam.h5", **criteria(filters))
    assert outputs_agree(tmp_path / "out.h5", tmp_path / "bam.h5")


def test_cram_decodes_against_the_reference_it_was_given(data, falsifiable, tmp_path):
    """Renames the reference recorded in the CRAM header, leaving cmuts hmm no
    source for the bases but --fasta."""
    converted = convert_format(data, tmp_path, CRAM)

    moved = tmp_path / "elsewhere.fasta"
    moved.write_bytes(converted.fasta.read_bytes())
    hidden = replace(converted, fasta=moved)

    run_cmuts(data, tmp_path / "bam.h5")

    with moved_aside(converted.fasta):
        summary = read_summary(run_cmuts(hidden, tmp_path / "cram.h5"))

    # The read counts agree whichever reference cmuts hmm decoded against, so
    # only the per-base fields depend on which one it used, and those are
    # written from the reads.
    falsifiable(summary.kept > 0)

    assert outputs_agree(tmp_path / "cram.h5", tmp_path / "bam.h5")
