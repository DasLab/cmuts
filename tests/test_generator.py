"""The generator itself, so that a bad fixture is reported as a bad fixture
rather than as a failing filter."""

import subprocess

import pytest

from datasets import DATASETS
from support import (
    generate, md_and_nm_tags, recomputed_md_and_nm_tags, records,
)


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_md_and_nm_match_samtools(datasets, tmp_path, name):
    """samtools recomputes both from nothing but the alignment and the reference, so
    agreement checks the generator against something other than itself."""
    data = datasets(name)

    assert md_and_nm_tags(data.bam) == recomputed_md_and_nm_tags(data, tmp_path)


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_alignments_are_well_formed(datasets, name):
    subprocess.run(["samtools", "quickcheck", str(datasets(name).bam)], check=True)


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_coordinate_sorted_without_sorting(datasets, name):
    subprocess.run(["samtools", "index", str(datasets(name).bam)],
                   check=True, capture_output=True)


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_header_lengths_match_the_reference(datasets, name):
    data = datasets(name)
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


def _divergences(data):
    """Reads whose stored sequence is longer than the reference span they
    align to, and reads whose is shorter.

    An insertion or a soft-clipped end stores bases that meet no reference
    position; a deletion meets positions no stored base covers.
    """
    longer = shorter = 0

    for line in records(data.bam, "-F", "4"):
        fields = line.split("\t")
        stored = len(fields[9])
        span = _reference_span(fields[5])

        longer += stored > span
        shorter += stored < span

    return longer, shorter


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_read_can_store_more_than_the_span_it_aligns_to(datasets, checked, name):
    checked(_divergences(datasets(name))[0])


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_read_can_store_less_than_the_span_it_aligns_to(datasets, checked, name):
    checked(_divergences(datasets(name))[1])


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
