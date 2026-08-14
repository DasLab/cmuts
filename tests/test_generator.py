"""The generator itself, so that a bad fixture is reported as a bad fixture and
not as a failing filter."""

import hashlib

import pytest

from alignments import NATIVE, generate
from oracle import (
    MAPPED_FLAG,
    header_lines,
    md_and_nm_tags,
    record_lines,
    records,
    recomputed_md_and_nm_tags,
    reference_lengths,
    reference_span,
    sequences,
)
from programs import samtools


@pytest.fixture(params=[NATIVE])
def fmt(request):
    """cmuts-gen writes one format, so its output is read in the format it
    wrote. A conversion here would test samtools."""
    return request.param


def _sq_fields(data) -> list:
    """Returns the fields of each @SQ line, keyed by tag."""
    return [
        dict(field.split(":", 1) for field in line.split("\t")[1:])
        for line in header_lines(data.bam, "@SQ")
    ]


def _alignment_count(data) -> int:
    """Returns how many records the file holds, mapped or not."""
    return data.mapped + data.unmapped


def test_checksums_match_hashlib(data, falsifiable):
    """hashlib computes the digests here, so agreement checks cmuts-gen against
    something other than the code cmuts-hmm checks them with."""
    written = {fields["SN"]: fields["M5"] for fields in _sq_fields(data)}

    falsifiable(len(written) > 0)

    assert written == {
        reference: hashlib.md5(seq.upper().encode()).hexdigest()
        for reference, seq in sequences(data.fasta).items()
    }


def test_md_and_nm_match_samtools(data, falsifiable, tmp_path):
    """samtools recomputes both from the alignment and the reference alone, so
    agreement checks the generator against something other than itself."""
    tags = md_and_nm_tags(data.bam)

    falsifiable(len(tags) > 0)

    assert tags == recomputed_md_and_nm_tags(data, tmp_path)


def test_alignments_are_well_formed(data, falsifiable):
    falsifiable(_alignment_count(data) > 0)

    samtools("quickcheck", data.bam)


def test_coordinate_sorted_without_sorting(data, falsifiable):
    falsifiable(_alignment_count(data) > 0)

    samtools("index", data.bam)


def test_header_lengths_match_the_reference(data, falsifiable):
    declared = {fields["SN"]: int(fields["LN"]) for fields in _sq_fields(data)}

    falsifiable(len(declared) > 0)

    assert declared == reference_lengths(data.fasta)


def test_the_same_seed_gives_the_same_data(tmp_path):
    first = generate(tmp_path, "a", seed=555, references=12, reads_per_ref=8)
    second = generate(tmp_path, "b", seed=555, references=12, reads_per_ref=8)

    assert first.fasta.read_bytes() == second.fasta.read_bytes()
    assert record_lines(first.bam) == record_lines(second.bam)


def test_a_different_seed_gives_different_data(tmp_path):
    first = generate(tmp_path, "c", seed=555, references=12, reads_per_ref=8)
    second = generate(tmp_path, "d", seed=556, references=12, reads_per_ref=8)

    assert record_lines(first.bam) != record_lines(second.bam)


def _reads_longer_than_their_span(data) -> int:
    """Returns how many reads store more bases than the reference span they
    align to, which an insertion or a soft-clipped end produces."""
    return sum(1 for record in records(data.bam, *MAPPED_FLAG)
               if len(record.sequence) > reference_span(record.cigar))


def _reads_shorter_than_their_span(data) -> int:
    """Returns how many reads store fewer bases than the reference span they
    align to, which a deletion produces."""
    return sum(1 for record in records(data.bam, *MAPPED_FLAG)
               if len(record.sequence) < reference_span(record.cigar))


def test_a_read_can_store_more_than_the_span_it_aligns_to(data, falsifiable):
    falsifiable(_reads_longer_than_their_span(data) > 0)


def test_a_read_can_store_less_than_the_span_it_aligns_to(data, falsifiable):
    falsifiable(_reads_shorter_than_their_span(data) > 0)
