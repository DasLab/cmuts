"""The generator itself, so that a bad fixture is reported as a bad fixture
rather than as a failing filter."""

import subprocess

import pytest

from conftest import SHAPES
from support import (
    Dataset, generate, md_and_nm_tags, recomputed_md_and_nm_tags, records,
)


@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_md_and_nm_match_samtools(datasets, tmp_path, shape):
    """samtools recomputes both from the alignment and the reference alone, so
    agreement checks the generator against something other than itself."""
    data = datasets(shape)

    assert md_and_nm_tags(data.bam) == recomputed_md_and_nm_tags(data, tmp_path)


@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_alignments_are_well_formed(datasets, shape):
    subprocess.run(["samtools", "quickcheck", str(datasets(shape).bam)], check=True)


@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_coordinate_sorted_without_sorting(datasets, shape):
    subprocess.run(["samtools", "index", str(datasets(shape).bam)],
                   check=True, capture_output=True)


@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_header_lengths_match_the_reference(datasets, shape):
    data = datasets(shape)
    declared = {}

    for line in records(data.bam, "-H"):
        if not line.startswith("@SQ"):
            continue
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:])
        declared[fields["SN"]] = int(fields["LN"])

    actual, name = {}, None
    for line in data.fasta.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            actual[name] = 0
        else:
            actual[name] += len(line)

    assert declared == actual


def test_the_same_seed_gives_the_same_data(tmp_path):
    first = generate(tmp_path, "a", seed=555, references=12, reads_per_ref=8)
    second = generate(tmp_path, "b", seed=555, references=12, reads_per_ref=8)

    assert first.fasta.read_bytes() == second.fasta.read_bytes()
    assert records(first.bam) == records(second.bam)


def test_a_different_seed_gives_different_data(tmp_path):
    first = generate(tmp_path, "c", seed=555, references=12, reads_per_ref=8)
    second = generate(tmp_path, "d", seed=556, references=12, reads_per_ref=8)

    assert records(first.bam) != records(second.bam)


def test_stored_length_diverges_from_aligned_span(datasets):
    longer = shorter = 0

    for line in records(datasets("plain").bam, "-F", "4"):
        fields = line.split("\t")
        stored = len(fields[9])
        span = _reference_span(fields[5])

        longer += stored > span
        shorter += stored < span

    assert longer > 0 and shorter > 0


def _reference_span(cigar):
    span, digits = 0, ""

    for character in cigar:
        if character.isdigit():
            digits += character
        else:
            if character in "MDN=X":
                span += int(digits)
            digits = ""

    return span
