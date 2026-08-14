"""The same alignments in each format htslib reads.

The format of a file decides how the bases and the CIGAR are stored and nothing
else, so a run over one format must answer as a run over another. CRAM stores a
sequence as its differences from a reference and so is the one format whose
reader needs a reference of its own; where it finds that reference is part of
the contract and is asserted here.
"""

import contextlib
import os
from dataclasses import replace

import pytest

from alignments import convert_format, replace_bases
from filters import COMPOUND, criteria, describe_filters
from oracle import assert_counts_agree
from outputs import outputs_agree
from programs import run_cmuts, try_cmuts

# The format the datasets are generated in, which the others are compared
# against.
NATIVE = "bam"

FORMATS = [NATIVE, "sam", "cram"]

# The criteria a format is tried under. The whole matrix belongs to
# test_filtering.py, which runs it over every dataset in the native format. A
# format can only change what a filter reads, so what is needed here is one run
# reading everything a filter can read and one run reading nothing.
CRITERIA = [{}, COMPOUND]


@contextlib.contextmanager
def moved_aside(path):
    """Renames a file for the duration of the block, so that nothing can read
    it."""
    aside = path.with_suffix(path.suffix + ".aside")
    os.rename(path, aside)
    try:
        yield
    finally:
        os.rename(aside, path)


@pytest.mark.parametrize("filters", CRITERIA, ids=describe_filters)
@pytest.mark.parametrize("fmt", FORMATS)
def test_every_format_gives_the_same_answer(data, falsifiable, tmp_path, fmt, filters):
    converted = convert_format(data, tmp_path, fmt)
    summary = run_cmuts(converted, tmp_path / "out.h5", **criteria(filters))

    # Only a run that kept reads compares what the format stored.
    falsifiable(summary.kept > 0)

    assert_counts_agree(summary, converted, criteria(filters))

    run_cmuts(data, tmp_path / "bam.h5", **criteria(filters))
    assert outputs_agree(tmp_path / "out.h5", tmp_path / "bam.h5")


def test_cram_decodes_against_the_reference_it_was_given(data, falsifiable, tmp_path):
    """Renames the reference recorded in the CRAM header, leaving cmuts-hmm no
    source for the bases but --fasta."""
    converted = convert_format(data, tmp_path, "cram")

    moved = tmp_path / "elsewhere.fasta"
    moved.write_bytes(converted.fasta.read_bytes())
    hidden = replace(converted, fasta=moved)

    run_cmuts(data, tmp_path / "bam.h5")

    with moved_aside(converted.fasta):
        summary = run_cmuts(hidden, tmp_path / "cram.h5")

    # The read counts agree whichever reference cmuts-hmm decoded against, so
    # only the per-base fields depend on which one it used, and those are
    # written from the reads.
    falsifiable(summary.kept > 0)

    assert outputs_agree(tmp_path / "cram.h5", tmp_path / "bam.h5")


def test_a_cram_is_refused_against_bases_it_was_not_written_from(data, falsifiable,
                                                                 tmp_path):
    """A CRAM stores its sequences as differences from a reference, so bases
    other than the ones it was written from fail to decode inside htslib. The
    same substitution on a BAM is caught by the M5 check instead, which
    test_reference.py covers.
    """
    converted = convert_format(data, tmp_path, "cram")

    # htslib fails while decoding a record, so a file holding no record never
    # reaches the decoder.
    falsifiable(data.mapped > 0)

    attempt = try_cmuts(replace_bases(converted, tmp_path), tmp_path / "out.h5")

    assert attempt.returncode != 0

    if data.mapped:
        assert "error reading alignment record" in attempt.stderr
