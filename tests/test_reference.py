"""The FASTA must hold the sequences the alignments were made against.

A name and a length describe a reference without identifying it: another
sequence of the same name and length passes every check that looks only at
shape, and its reads are then scored against the wrong bases. An @SQ M5
checksum is taken over the bases themselves, for the references whose aligner
wrote one.
"""

import hashlib
from dataclasses import replace

import pytest

from datasets import DATASETS
from support import (
    outputs_agree, reheadered, run_cmuts, sequences, substituted, try_cmuts, written,
)


def has_references(data) -> bool:
    """Whether the FASTA declares a reference to take a checksum over."""
    return len(sequences(data.fasta)) > 0


def sampled(data):
    """The first reference the FASTA declares, the last, and one between."""
    names = sorted(sequences(data.fasta))

    return [names[0], names[len(names) // 2], names[-1]] if names else []


# ---------------------------------------------------------------------------
# Rewriting a reference and the checksum its header declares
# ---------------------------------------------------------------------------


def md5(seq):
    """The digest of a sequence, computed by something other than the code
    under test."""
    return hashlib.md5(seq.upper().encode()).hexdigest()


def with_checksums(data, tmp_path, checksum=md5, only=None):
    """The same alignments, their header declaring an M5 for each reference."""
    records = sequences(data.fasta)

    def declare(line):
        name = next(field[3:] for field in line.split("\t") if field.startswith("SN:"))

        if only is not None and name not in only:
            return line

        return line + "\tM5:" + checksum(records[name])

    def rewrite(text):
        return "".join(
            (declare(line) if line.startswith("@SQ") else line) + "\n"
            for line in text.splitlines()
        )

    return reheadered(data, tmp_path, rewrite)


# ---------------------------------------------------------------------------
# What a checksum settles
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_matching_checksum_is_accepted(datasets, falsifiable, tmp_path, name):
    data = datasets(name)

    falsifiable(has_references(data))

    checked = run_cmuts(with_checksums(data, tmp_path), tmp_path / "checked.h5")
    plain = run_cmuts(data, tmp_path / "plain.h5")

    assert checked == plain
    assert outputs_agree(tmp_path / "checked.h5", tmp_path / "plain.h5")


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_substituted_reference_is_refused(datasets, falsifiable, tmp_path, name):
    data = datasets(name)

    falsifiable(has_references(data))

    wrong = substituted(with_checksums(data, tmp_path), tmp_path)
    attempt = try_cmuts(wrong, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not the sequence the alignments were made against" in attempt.stderr


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_reference_without_a_checksum_is_not_checked(datasets, falsifiable, tmp_path,
                                                       name):
    """Runs to the end against bases the alignments were not made from, and
    scores them: the result a checksum exists to refuse."""
    data = datasets(name)
    wrong, right = tmp_path / "wrong.h5", tmp_path / "right.h5"

    assert try_cmuts(substituted(data, tmp_path), wrong).returncode == 0

    kept = run_cmuts(data, right).kept

    # Where every read was rejected the two runs count nothing either way, so
    # the bases they counted it over cannot be told apart.
    falsifiable(kept > 0)

    if kept:
        assert not outputs_agree(wrong, right), \
            "the substituted bases are not the ones that were scored"


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_the_checksum_checked_is_the_one_for_that_reference(datasets, falsifiable,
                                                            tmp_path, name):
    data = datasets(name)
    declared = with_checksums(data, tmp_path)
    reached = sampled(data)

    falsifiable(len(reached) > 0)

    for reference in reached:
        attempt = try_cmuts(substituted(declared, tmp_path, only={reference}),
                            tmp_path / "out.h5", overwrite=True)

        assert attempt.returncode != 0, reference
        assert f'"{reference}"' in attempt.stderr


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_checksum_on_a_single_reference_is_still_checked(datasets, falsifiable,
                                                           tmp_path, name):
    data = datasets(name)
    last = sampled(data)[-1:]

    falsifiable(len(last) > 0)

    for reference in last:
        declared = with_checksums(data, tmp_path, only={reference})
        attempt = try_cmuts(substituted(declared, tmp_path, only={reference}),
                            tmp_path / "out.h5")

        assert attempt.returncode != 0
        assert f'"{reference}"' in attempt.stderr


# ---------------------------------------------------------------------------
# Reading the field it is written in
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_checksum_is_read_whatever_its_case(datasets, falsifiable, tmp_path, name):
    """The SAM spec fixes M5 as hexadecimal without fixing its case."""
    data = datasets(name)

    falsifiable(has_references(data))

    shouting = with_checksums(data, tmp_path, checksum=lambda seq: md5(seq).upper())

    assert try_cmuts(shouting, tmp_path / "out.h5").returncode == 0


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_soft_masked_reference_hashes_as_an_upper_case_one(datasets, falsifiable,
                                                             tmp_path, name):
    """The SAM spec defines M5 over the sequence uppercased."""
    data = datasets(name)

    falsifiable(has_references(data))

    masked = replace(data, fasta=written(
        {name: seq.lower() for name, seq in sequences(data.fasta).items()},
        tmp_path / "masked.fasta"))

    assert try_cmuts(with_checksums(masked, tmp_path), tmp_path / "out.h5").returncode == 0


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_checksum_is_read_only_to_the_end_of_its_field(datasets, falsifiable,
                                                         tmp_path, name):
    data = datasets(name)

    falsifiable(has_references(data))

    trailing = with_checksums(data, tmp_path,
                              checksum=lambda seq: md5(seq) + "\tUR:file:/nowhere")

    assert try_cmuts(trailing, tmp_path / "out.h5").returncode == 0


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_checksum_that_is_not_an_md5_is_refused(datasets, falsifiable, tmp_path, name):
    data = datasets(name)

    falsifiable(has_references(data))

    truncated = with_checksums(data, tmp_path, checksum=lambda seq: md5(seq)[:8])
    attempt = try_cmuts(truncated, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not an MD5" in attempt.stderr
