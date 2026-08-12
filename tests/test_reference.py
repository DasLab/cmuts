"""The FASTA must hold the sequences the alignments were made against.

A name and a length describe a reference without identifying it: another
sequence of the same name and length passes every check on shape alone, and
its reads are then scored against the wrong bases. An @SQ M5 checksum is taken
over the bases themselves, for the references whose aligner wrote one.
"""

import hashlib
from dataclasses import replace

import pytest

from support import outputs_agree, reheadered, run_cmuts, sequences, try_cmuts

COMPLEMENT = str.maketrans("ACGT", "TGCA")


# ---------------------------------------------------------------------------
# Rewriting a reference, and what the header says about it
# ---------------------------------------------------------------------------


def written(records, path):
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()))
    return path


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


def substituted(data, tmp_path, only=None):
    """A FASTA agreeing on every name and length, holding different bases."""
    records = {
        name: seq.translate(COMPLEMENT) if only is None or name in only else seq
        for name, seq in sequences(data.fasta).items()
    }

    return replace(data, fasta=written(records, tmp_path / "substituted.fasta"))


# ---------------------------------------------------------------------------
# What a checksum settles
# ---------------------------------------------------------------------------


def test_a_matching_checksum_is_accepted(data, tmp_path):
    """Verification decides whether a run happens and not what it produces, so
    the result matches a run with no checksums declared."""
    checked = run_cmuts(with_checksums(data, tmp_path), tmp_path / "checked.h5")
    plain = run_cmuts(data, tmp_path / "plain.h5")

    assert checked == plain
    assert outputs_agree(tmp_path / "checked.h5", tmp_path / "plain.h5")


def test_a_substituted_reference_is_refused(data, tmp_path):
    """Same names, same lengths, different bases."""
    wrong = substituted(with_checksums(data, tmp_path), tmp_path)
    attempt = try_cmuts(wrong, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not the sequence the alignments were made against" in attempt.stderr


def test_a_reference_without_a_checksum_is_not_checked(data, tmp_path):
    """With no M5 there is nothing to compare against, so a substituted FASTA
    of the right shape still runs."""
    assert try_cmuts(substituted(data, tmp_path), tmp_path / "out.h5").returncode == 0


def test_the_checksum_checked_is_the_one_for_that_reference(data, tmp_path):
    """Every reference declares a checksum and exactly one sequence is
    replaced, so a comparison a line out of step accepts a substituted
    reference or refuses an untouched one."""
    declared = with_checksums(data, tmp_path)
    names = list(sequences(data.fasta))

    for name in (names[0], names[len(names) // 2], names[-1]):
        attempt = try_cmuts(substituted(declared, tmp_path, only={name}),
                            tmp_path / "out.h5", overwrite=True)

        assert attempt.returncode != 0, name
        assert f'"{name}"' in attempt.stderr


def test_a_checksum_on_one_reference_alone_is_still_checked(data, tmp_path):
    """Most @SQ lines declare no M5, so the one that does must still be found
    among them."""
    last = list(sequences(data.fasta))[-1]
    declared = with_checksums(data, tmp_path, only={last})

    attempt = try_cmuts(substituted(declared, tmp_path, only={last}), tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert f'"{last}"' in attempt.stderr


# ---------------------------------------------------------------------------
# Reading the field it is written in
# ---------------------------------------------------------------------------


def test_a_checksum_is_read_whatever_its_case(data, tmp_path):
    """The spec fixes hexadecimal, not which case it is spelt in."""
    shouting = with_checksums(data, tmp_path, checksum=lambda seq: md5(seq).upper())

    assert try_cmuts(shouting, tmp_path / "out.h5").returncode == 0


def test_a_soft_masked_reference_hashes_as_an_upper_case_one(data, tmp_path):
    """Lower case in a FASTA marks repeats rather than different bases, and the
    digest is defined over the sequence uppercased, so masking must not change
    whether a run is refused."""
    masked = replace(data, fasta=written(
        {name: seq.lower() for name, seq in sequences(data.fasta).items()},
        tmp_path / "masked.fasta"))

    assert try_cmuts(with_checksums(masked, tmp_path), tmp_path / "out.h5").returncode == 0


def test_a_checksum_is_read_only_to_the_end_of_its_field(data, tmp_path):
    """M5 need not be the last tag on its line: the value ends where its field
    does, not where the line does."""
    trailing = with_checksums(data, tmp_path,
                              checksum=lambda seq: md5(seq) + "\tUR:file:/nowhere")

    assert try_cmuts(trailing, tmp_path / "out.h5").returncode == 0


def test_a_checksum_that_is_not_an_md5_is_refused(data, tmp_path):
    """A value that cannot be a digest is refused rather than skipped, which
    would report a check that was never made."""
    truncated = with_checksums(data, tmp_path, checksum=lambda seq: md5(seq)[:8])
    attempt = try_cmuts(truncated, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not an MD5" in attempt.stderr
