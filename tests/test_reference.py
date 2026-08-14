"""The FASTA must hold the sequences the alignments were made against.

A name and a length describe a reference without identifying it: another
sequence of the same name and length passes every check that looks only at
shape, and its reads are then scored against the wrong bases. An @SQ M5
checksum is taken over the bases themselves, for the references whose aligner
wrote one.
"""

from dataclasses import replace

import pytest

from datasets import DATASETS
from support import (
    outputs_agree, rename_references, replace_bases, replace_checksums, run_cmuts,
    sequences, try_cmuts, write_fasta,
)


def has_references(data) -> bool:
    """Whether the FASTA holds a reference to take a checksum over."""
    return len(sequences(data.fasta)) > 0


def sample_references(data):
    """The first reference the FASTA holds, the last, and one between."""
    names = sorted(sequences(data.fasta))

    return [names[0], names[len(names) // 2], names[-1]] if names else []


# ---------------------------------------------------------------------------
# What a checksum settles
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_matching_checksum_changes_nothing_counted(datasets, falsifiable, tmp_path,
                                                     name):
    """cmuts-gen writes a matching checksum into every header, so every run in
    the suite already asserts that one is accepted. What is left to state is
    that checking it decides nothing about the result."""
    data = datasets(name)

    falsifiable(has_references(data))

    checked = run_cmuts(data, tmp_path / "checked.h5")
    unchecked = run_cmuts(data, tmp_path / "unchecked.h5", verify="name,length")

    assert checked == unchecked
    assert outputs_agree(tmp_path / "checked.h5", tmp_path / "unchecked.h5")


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_substituted_reference_is_refused(datasets, falsifiable, tmp_path, name):
    data = datasets(name)

    falsifiable(has_references(data))

    attempt = try_cmuts(replace_bases(data, tmp_path), tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not the sequence the alignments were made against" in attempt.stderr


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_reference_without_a_checksum_is_not_checked(datasets, falsifiable, tmp_path,
                                                       name):
    """Runs to the end against bases the alignments were not made from, and
    scores them: the result a checksum exists to refuse."""
    data = datasets(name)
    wrong, right = tmp_path / "wrong.h5", tmp_path / "right.h5"
    stripped = replace_checksums(data, tmp_path, lambda reference, m5: None)

    assert try_cmuts(replace_bases(stripped, tmp_path), wrong).returncode == 0

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
    reached = sample_references(data)

    falsifiable(len(reached) > 0)

    for reference in reached:
        attempt = try_cmuts(replace_bases(data, tmp_path, only={reference}),
                            tmp_path / "out.h5", overwrite=True)

        assert attempt.returncode != 0, reference
        assert f'"{reference}"' in attempt.stderr


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_checksum_on_a_single_reference_is_still_checked(datasets, falsifiable,
                                                           tmp_path, name):
    data = datasets(name)
    last = sample_references(data)[-1:]

    falsifiable(len(last) > 0)

    for kept in last:
        alone = replace_checksums(
            data, tmp_path, lambda reference, m5: m5 if reference == kept else None)
        attempt = try_cmuts(replace_bases(alone, tmp_path, only={kept}),
                            tmp_path / "out.h5")

        assert attempt.returncode != 0
        assert f'"{kept}"' in attempt.stderr


# ---------------------------------------------------------------------------
# Choosing which comparisons to make
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_names_go_unchecked_where_they_are_left_out(datasets, falsifiable, tmp_path,
                                                    name):
    """A FASTA renamed after the alignments were made is the reason to leave the
    name out, so the counts must come to what the original names give."""
    data = datasets(name)

    falsifiable(has_references(data))

    against_names = rename_references(data, tmp_path)
    attempt = try_cmuts(against_names, tmp_path / "renamed.h5",
                        verify="length,checksum")

    assert attempt.returncode == 0, attempt.stderr

    run_cmuts(data, tmp_path / "plain.h5")

    assert outputs_agree(tmp_path / "renamed.h5", tmp_path / "plain.h5")


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_checksum_goes_unchecked_where_it_is_left_out(datasets, falsifiable,
                                                        tmp_path, name):
    data = datasets(name)

    falsifiable(has_references(data))

    wrong = replace_bases(data, tmp_path)

    assert try_cmuts(wrong, tmp_path / "out.h5", verify="name,length").returncode == 0


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_verifying_none_takes_the_fasta_on_trust(datasets, falsifiable, tmp_path, name):
    """The FASTA disagrees on every count a check could make: cmuts-gen wrote
    the checksums over the bases it held before either was changed."""
    data = datasets(name)

    falsifiable(has_references(data))

    wrong = replace_bases(rename_references(data, tmp_path), tmp_path)

    assert try_cmuts(wrong, tmp_path / "out.h5", verify="none").returncode == 0


# ---------------------------------------------------------------------------
# Reading the field it is written in
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_checksum_is_read_whatever_its_case(datasets, falsifiable, tmp_path, name):
    """The SAM spec fixes M5 as hexadecimal without fixing its case."""
    data = datasets(name)

    falsifiable(has_references(data))

    shouting = replace_checksums(data, tmp_path, lambda reference, m5: m5.upper())

    assert try_cmuts(shouting, tmp_path / "out.h5").returncode == 0


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_soft_masked_reference_hashes_as_an_upper_case_one(datasets, falsifiable,
                                                             tmp_path, name):
    """The SAM spec defines M5 over the sequence uppercased."""
    data = datasets(name)

    falsifiable(has_references(data))

    masked = replace(data, fasta=write_fasta(
        {reference: seq.lower() for reference, seq in sequences(data.fasta).items()},
        tmp_path / "masked.fasta"))

    assert try_cmuts(masked, tmp_path / "out.h5").returncode == 0


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_checksum_is_read_only_to_the_end_of_its_field(datasets, falsifiable,
                                                         tmp_path, name):
    data = datasets(name)

    falsifiable(has_references(data))

    trailing = replace_checksums(data, tmp_path,
                                 lambda reference, m5: m5 + "\tUR:file:/nowhere")

    assert try_cmuts(trailing, tmp_path / "out.h5").returncode == 0


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_checksum_that_is_not_an_md5_is_refused(datasets, falsifiable, tmp_path, name):
    data = datasets(name)

    falsifiable(has_references(data))

    truncated = replace_checksums(data, tmp_path, lambda reference, m5: m5[:8])
    attempt = try_cmuts(truncated, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not an MD5" in attempt.stderr
