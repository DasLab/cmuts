"""The FASTA must hold the sequences the alignments were made against.

A name and a length describe a reference without identifying it: another
sequence of the same name and length passes every check that looks only at
shape, and its reads are then scored against the wrong bases. An @SQ M5
checksum is computed over the bases themselves.
"""

from dataclasses import replace

import pytest

from alignments import (
    SELF_CONTAINED,
    lengthen_references,
    rename_references,
    replace_bases,
    replace_checksums,
    write_fasta,
)
from oracle import sequences
from outputs import outputs_agree, read_summary
from programs import run_cmuts, try_cmuts


def has_references(data) -> bool:
    """Returns whether the FASTA holds a reference to take a checksum over."""
    return len(sequences(data.fasta)) > 0


def sample_references(data) -> list:
    """Returns the first reference of the FASTA, the last, and one in
    between."""
    names = sorted(sequences(data.fasta))

    return [names[0], names[len(names) // 2], names[-1]] if names else []


# ---------------------------------------------------------------------------
# What a checksum settles
# ---------------------------------------------------------------------------


def test_a_matching_checksum_changes_nothing_counted(data, falsifiable, tmp_path):
    """cmuts gen writes a matching checksum into every header, so every run in
    the suite already asserts that a matching one is accepted. This test
    asserts that checking it leaves the result unchanged."""
    falsifiable(has_references(data))

    checked = read_summary(run_cmuts(data, tmp_path / "checked.h5"))
    unchecked = read_summary(
        run_cmuts(data, tmp_path / "unchecked.h5", verify="name"))

    assert checked == unchecked
    assert outputs_agree(tmp_path / "checked.h5", tmp_path / "unchecked.h5")


def test_a_substituted_reference_is_refused(data, falsifiable, tmp_path):
    falsifiable(has_references(data))

    attempt = try_cmuts(replace_bases(data, tmp_path), tmp_path / "out.h5")

    assert attempt.returncode != 0


@pytest.mark.parametrize("fmt", SELF_CONTAINED, indirect=True)
def test_a_reference_without_a_checksum_is_not_checked(data, falsifiable, tmp_path):
    """Runs to the end against bases the alignments were not made from, and
    scores them: the result a checksum exists to prevent.

    Reading the wrong bases at all requires a format that stores its own
    sequence.
    """
    wrong, right = tmp_path / "wrong.h5", tmp_path / "right.h5"
    stripped = replace_checksums(data, tmp_path, lambda reference, m5: None)

    assert try_cmuts(replace_bases(stripped, tmp_path), wrong).returncode == 0

    kept = read_summary(run_cmuts(data, right)).kept

    # Where every read was rejected, the two runs count nothing either way, so
    # the bases they would have counted over cannot be told apart.
    falsifiable(kept > 0)

    if kept:
        assert not outputs_agree(wrong, right), \
            "the substituted bases are not the ones that were scored"


def test_a_substitution_in_any_reference_is_refused(data, falsifiable, tmp_path):
    """Substitutes the first reference of the FASTA, the last, and one in
    between, so that a check reaching only the first is caught."""
    reached = sample_references(data)

    falsifiable(len(reached) > 0)

    for reference in reached:
        attempt = try_cmuts(replace_bases(data, tmp_path, only={reference}),
                            tmp_path / "out.h5", overwrite=True)

        assert attempt.returncode != 0, reference


def test_a_checksum_on_a_single_reference_is_still_checked(data, falsifiable, tmp_path):
    last = sample_references(data)[-1:]

    falsifiable(len(last) > 0)

    for kept in last:
        alone = replace_checksums(
            data, tmp_path, lambda reference, m5: m5 if reference == kept else None)
        attempt = try_cmuts(replace_bases(alone, tmp_path, only={kept}),
                            tmp_path / "out.h5")

        assert attempt.returncode != 0


# ---------------------------------------------------------------------------
# Choosing which comparisons to make
# ---------------------------------------------------------------------------


def test_names_go_unchecked_where_they_are_left_out(data, falsifiable, tmp_path):
    """A FASTA renamed after the alignments were made is the reason to leave
    the name out, so the counts must come to what the original names give."""
    falsifiable(has_references(data))

    against_names = rename_references(data, tmp_path)
    attempt = try_cmuts(against_names, tmp_path / "renamed.h5",
                        verify="checksum")

    assert attempt.returncode == 0, attempt.stderr

    run_cmuts(data, tmp_path / "plain.h5")

    assert outputs_agree(tmp_path / "renamed.h5", tmp_path / "plain.h5")


@pytest.mark.parametrize("fmt", SELF_CONTAINED, indirect=True)
def test_a_checksum_goes_unchecked_where_it_is_left_out(data, falsifiable, tmp_path):
    """Reading the wrong bases at all requires a format that stores its own
    sequence."""
    falsifiable(has_references(data))

    wrong = replace_bases(data, tmp_path)

    assert try_cmuts(wrong, tmp_path / "out.h5", verify="name").returncode == 0


def test_verifying_none_takes_the_fasta_on_trust(data, falsifiable, tmp_path):
    """The FASTA disagrees on every count a check could make: cmuts gen wrote
    the checksums over the bases it held before either was changed."""
    falsifiable(has_references(data))

    wrong = replace_bases(rename_references(data, tmp_path), tmp_path)

    assert try_cmuts(wrong, tmp_path / "out.h5", verify="none").returncode == 0


@pytest.mark.parametrize("verify", ["name,checksum", "name", "checksum", "none"])
def test_a_length_the_header_does_not_declare_is_refused(data, falsifiable, tmp_path,
                                                         verify):
    """A record must be the length its header declares, and --verify cannot turn
    that off. The per-reference buffers are sized from the declared lengths
    before any record is read, so a longer record would overrun one."""
    falsifiable(has_references(data))

    longer = lengthen_references(data, tmp_path)

    assert try_cmuts(longer, tmp_path / "out.h5", verify=verify).returncode != 0


# ---------------------------------------------------------------------------
# Reading the field the checksum is written in
# ---------------------------------------------------------------------------


def test_a_checksum_is_read_whatever_its_case(data, falsifiable, tmp_path):
    """The SAM spec fixes M5 as hexadecimal without fixing its case."""
    falsifiable(has_references(data))

    shouting = replace_checksums(data, tmp_path, lambda reference, m5: m5.upper())

    assert try_cmuts(shouting, tmp_path / "out.h5").returncode == 0


def test_a_soft_masked_reference_hashes_as_an_upper_case_one(data, falsifiable,
                                                             tmp_path):
    """The SAM spec defines M5 over the sequence uppercased."""
    falsifiable(has_references(data))

    masked = replace(data, fasta=write_fasta(
        {reference: seq.lower() for reference, seq in sequences(data.fasta).items()},
        tmp_path / "masked.fasta"))

    assert try_cmuts(masked, tmp_path / "out.h5").returncode == 0


def test_a_checksum_is_read_only_to_the_end_of_its_field(data, falsifiable, tmp_path):
    falsifiable(has_references(data))

    trailing = replace_checksums(data, tmp_path,
                                 lambda reference, m5: m5 + "\tUR:file:/nowhere")

    assert try_cmuts(trailing, tmp_path / "out.h5").returncode == 0


def test_a_checksum_that_is_not_an_md5_is_refused(data, falsifiable, tmp_path):
    falsifiable(has_references(data))

    truncated = replace_checksums(data, tmp_path, lambda reference, m5: m5[:8])
    attempt = try_cmuts(truncated, tmp_path / "out.h5")

    assert attempt.returncode != 0
