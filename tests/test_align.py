"""Aligning reads with cmuts align.

The reads are taken from an alignment that cmuts gen wrote and read back as
FASTQ, so every test has a reference that its reads align to.
"""

from __future__ import annotations

import shutil

import pytest

from alignments import NATIVE, measure_dataset
from oracle import sequences
from programs import (CMUTS_ALIGN, attempt, run_align, samtools, samtools_into,
                      try_align, try_cmuts)

# A BAM file stores the length of a read name in one byte.
NAME_LIMIT = 254

# The dataset the reads are taken from. One is enough, because no contract here
# depends on how the reads are laid out.
DATASET = "plain"

# The preset every test aligns under, other than those naming their own.
PRESET = "sr"

# How much of a read each mate holds when one read is written out as a pair.
# Anything above half leaves the two mates overlapping, which the merger
# requires.
MATE_SHARE = 0.7

pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None,
                                reason="minimap2 is required to align reads")


# ---------------------------------------------------------------------------
# Reading the output
# ---------------------------------------------------------------------------


def header_of(bam) -> str:
    return samtools("view", "-H", bam)


def reference_lines(bam) -> list:
    return [line for line in header_of(bam).splitlines() if line.startswith("@SQ")]


def read_names(bam) -> list:
    return [line.split("\t")[0] for line in samtools("view", bam).splitlines()]


def checksums_in(text: str) -> set:
    """The M5 field of every @SQ line, whether of a header or of a dictionary."""
    return {field for line in text.splitlines() if line.startswith("@SQ")
            for field in line.split("\t") if field.startswith("M5:")}


# ---------------------------------------------------------------------------
# Writing the reads
# ---------------------------------------------------------------------------


def rewrite_names(reads, destination, rename):
    """Copies a FASTQ file with every read name replaced."""
    lines = reads.read_text().splitlines()

    with open(destination, "w") as handle:
        for at in range(0, len(lines), 4):
            handle.write("@" + rename(lines[at][1:]) + "\n")
            handle.write("\n".join(lines[at + 1:at + 4]) + "\n")

    return destination


def unaligned_bam(reads, destination):
    """Writes the reads out as an unaligned BAM file."""
    samtools("import", "-0", reads, "-o", destination)

    return destination


def split_into_mates(reads, directory):
    """Writes a FASTQ file out as two, each mate holding most of one read, so
    that the two mates overlap."""
    complement = str.maketrans("ACGTN", "TGCAN")
    lines = reads.read_text().splitlines()
    first = directory / "r1.fastq"
    second = directory / "r2.fastq"

    with open(first, "w") as one, open(second, "w") as two:
        for at in range(0, len(lines), 4):
            name, bases, _, qualities = lines[at:at + 4]
            take = int(len(bases) * MATE_SHARE)

            one.write(f"{name}\n{bases[:take]}\n+\n{qualities[:take]}\n")
            two.write(f"{name}\n{bases[-take:][::-1].translate(complement)}\n+\n"
                      f"{qualities[-take:][::-1]}\n")

    return first, second


@pytest.fixture(scope="session")
def alignable(catalogue, tmp_path_factory):
    """A FASTA file, and a FASTQ file holding reads that align to it."""
    data = catalogue(DATASET, NATIVE)
    reads = tmp_path_factory.mktemp("reads") / "reads.fastq"

    samtools_into(reads, "fastq", data.bam)

    return data.fasta, reads


@pytest.fixture
def aligned(alignable, tmp_path):
    """The reads aligned to the reference, as one sorted BAM file."""
    fasta, reads = alignable

    return run_align(fasta, tmp_path / "aligned.bam", PRESET, reads)


# ---------------------------------------------------------------------------
# What the output holds
# ---------------------------------------------------------------------------


def test_the_output_is_coordinate_sorted(aligned):
    assert "SO:coordinate" in header_of(aligned)


def test_the_output_holds_the_reads_that_aligned(aligned):
    assert len(read_names(aligned)) > 0


def test_every_reference_declares_the_checksum_of_its_sequence(alignable, aligned):
    fasta, _ = alignable
    declared = checksums_in(header_of(aligned))

    assert len(declared) == len(reference_lines(aligned))
    assert declared <= checksums_in(samtools("dict", fasta))


def test_a_read_name_too_long_for_a_bam_file_is_shortened(alignable, tmp_path):
    fasta, reads = alignable
    renamed = rewrite_names(reads, tmp_path / "long.fastq",
                            lambda name: name + "x" * NAME_LIMIT)

    aligned = run_align(fasta, tmp_path / "long.bam", PRESET, renamed)

    assert read_names(aligned)
    assert max(len(name) for name in read_names(aligned)) == NAME_LIMIT


def test_the_output_can_be_counted(alignable, aligned, tmp_path):
    fasta, _ = alignable
    counted = measure_dataset("aligned", (aligned,), fasta)

    assert try_cmuts(counted, tmp_path / "counts.h5").returncode == 0


# ---------------------------------------------------------------------------
# The forms the reads arrive in
# ---------------------------------------------------------------------------


def test_the_same_reads_align_alike_from_a_bam_and_from_a_fastq(
        alignable, aligned, tmp_path):
    fasta, reads = alignable
    unaligned = unaligned_bam(reads, tmp_path / "unaligned.bam")

    from_bam = run_align(fasta, tmp_path / "frombam.bam", PRESET, unaligned)

    assert samtools("view", from_bam) == samtools("view", aligned)


@pytest.mark.skipif(shutil.which("vsearch") is None,
                    reason="vsearch is required to merge a pair of mates")
def test_a_pair_of_mates_is_merged_and_aligned(alignable, tmp_path):
    fasta, reads = alignable
    first, second = split_into_mates(reads, tmp_path)

    merged = run_align(fasta, tmp_path / "merged.bam", PRESET, first, second)

    assert len(read_names(merged)) > 0


# ---------------------------------------------------------------------------
# What is refused
# ---------------------------------------------------------------------------


def test_an_already_aligned_bam_is_refused(catalogue, tmp_path):
    data = catalogue(DATASET, NATIVE)

    assert try_align(data.fasta, tmp_path / "out.bam", PRESET, data.bam).returncode != 0


def test_a_bam_given_as_one_file_of_a_pair_is_refused(alignable, tmp_path):
    fasta, reads = alignable
    unaligned = unaligned_bam(reads, tmp_path / "unaligned.bam")

    assert try_align(fasta, tmp_path / "out.bam", PRESET,
                     reads, unaligned).returncode != 0


def test_a_pair_of_mates_is_refused_under_a_long_read_preset(alignable, tmp_path):
    fasta, reads = alignable
    first, second = split_into_mates(reads, tmp_path)

    assert try_align(fasta, tmp_path / "out.bam", "map-ont",
                     first, second).returncode != 0


def test_a_preset_that_aligns_reads_against_each_other_is_refused(alignable, tmp_path):
    fasta, reads = alignable

    assert try_align(fasta, tmp_path / "out.bam", "ava-ont", reads).returncode != 0


def test_a_truncated_bam_is_refused(alignable, tmp_path):
    fasta, reads = alignable
    whole = unaligned_bam(reads, tmp_path / "unaligned.bam").read_bytes()
    truncated = tmp_path / "truncated.bam"
    output = tmp_path / "out.bam"

    truncated.write_bytes(whole[:len(whole) // 2])

    assert try_align(fasta, output, PRESET, truncated).returncode != 0
    assert not output.exists()


def test_a_file_holding_neither_reads_nor_alignments_is_refused(alignable, tmp_path):
    fasta, _ = alignable

    assert try_align(fasta, tmp_path / "out.bam", PRESET, fasta).returncode != 0


def test_a_minimap2_index_in_place_of_the_fasta_is_refused(alignable, tmp_path):
    fasta, reads = alignable
    index = tmp_path / "ref.mmi"

    attempt(["minimap2", "-d", index, fasta])

    assert try_align(index, tmp_path / "out.bam", PRESET, reads).returncode != 0


def test_a_fasta_that_repeats_a_reference_name_is_refused(alignable, tmp_path):
    """The report names the file and the reference that repeats, rather than
    the failure that the repeated name causes further down the pipeline."""
    fasta, reads = alignable
    repeated = tmp_path / "repeated.fasta"
    name = next(iter(sequences(fasta)))

    repeated.write_text(fasta.read_text() * 2)
    result = try_align(repeated, tmp_path / "out.bam", PRESET, reads)
    reason = result.stderr.strip().splitlines()[-1]

    assert result.returncode != 0
    assert str(repeated) in reason
    assert f'"{name}"' in reason


def test_an_existing_output_is_kept_unless_overwriting_is_asked_for(alignable, aligned):
    fasta, reads = alignable
    before = aligned.read_bytes()

    assert try_align(fasta, aligned, PRESET, reads).returncode != 0
    assert aligned.read_bytes() == before

    run_align(fasta, aligned, PRESET, reads, overwrite=True)


def test_a_directory_in_place_of_the_output_is_refused(alignable, tmp_path):
    fasta, reads = alignable

    assert try_align(fasta, tmp_path, PRESET, reads).returncode != 0


def test_a_run_that_is_refused_leaves_no_partial_file(alignable, tmp_path):
    fasta, _ = alignable
    output = tmp_path / "out.bam"

    assert try_align(fasta, output, PRESET, fasta).returncode != 0
    assert list(tmp_path.glob("out.bam*")) == []


def test_the_program_reports_its_version():
    assert attempt([*CMUTS_ALIGN, "--version"]).returncode == 0
