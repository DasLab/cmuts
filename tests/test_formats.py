"""The same alignments in each format htslib reads.

Which format a file is in decides how the bases and the CIGAR are stored and
nothing else, so a run over one must answer as a run over another. CRAM stores
the sequence as differences from a reference and so is the one format whose
reader needs a reference of its own; where it finds one is part of the contract
and is asserted here.
"""

import contextlib
import os
from dataclasses import replace

import pytest

from datasets import DATASETS
from filters import COMPOUND, criteria, describe_filters
from support import (
    assert_counts_agree, converted, outputs_agree, run_cmuts, substituted, try_cmuts,
)

# The format the datasets are generated in, which the others are compared
# against.
NATIVE = "bam"

FORMATS = [NATIVE, "sam", "cram"]

# Which criteria a format is tried under. The whole matrix belongs to
# test_filtering.py, which runs it over every dataset in the native format; a
# format can only change what a filter reads, so what is wanted here is a run
# reading everything a filter can read and a run reading nothing.
CRITERIA = [{}, COMPOUND]


@contextlib.contextmanager
def unreachable(path):
    """Renames a file aside for the duration, so nothing can fall back to it."""
    aside = path.with_suffix(path.suffix + ".aside")
    os.rename(path, aside)
    try:
        yield
    finally:
        os.rename(aside, path)


@pytest.mark.parametrize("filters", CRITERIA, ids=describe_filters)
@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("name", sorted(DATASETS))
def test_every_format_gives_the_same_answer(datasets, tmp_path, name, fmt, filters):
    native = datasets(name)
    data = converted(native, tmp_path, fmt)
    summary = run_cmuts(data, tmp_path / "out.h5", **criteria(filters))

    assert_counts_agree(summary, data, criteria(filters))

    run_cmuts(native, tmp_path / "bam.h5", **criteria(filters))
    assert outputs_agree(tmp_path / "out.h5", tmp_path / "bam.h5")


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_cram_decodes_against_the_reference_it_was_given(datasets, tmp_path, name):
    """Renames the reference recorded in the CRAM header aside, leaving cmuts-hmm
    no source for the bases but --fasta."""
    native = datasets(name)
    data = converted(native, tmp_path, "cram")

    moved = tmp_path / "elsewhere.fasta"
    moved.write_bytes(data.fasta.read_bytes())
    hidden = replace(data, fasta=moved)

    run_cmuts(native, tmp_path / "bam.h5")

    with unreachable(data.fasta):
        run_cmuts(hidden, tmp_path / "cram.h5")

    # Read counts agree whichever reference cmuts-hmm decoded against, so only the
    # per-base fields show which one it used.
    assert outputs_agree(tmp_path / "cram.h5", tmp_path / "bam.h5")


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_cram_is_refused_against_bases_it_was_not_written_from(datasets, tmp_path,
                                                                 name):
    """A CRAM stores its sequences as differences from a reference and carries
    an M5 for each, samtools declaring one as it converts. Decoding it against
    other bases fails in htslib, before the checksum in the header is reached,
    so this asserts the refusal and not the reason given for it.
    """
    data = converted(datasets(name), tmp_path, "cram")

    attempt = try_cmuts(substituted(data, tmp_path), tmp_path / "out.h5")

    assert attempt.returncode != 0
