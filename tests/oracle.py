"""Tools for deriving the expected result of a run from its inputs.

Every value here comes from the alignment and the FASTA alone, by way of
samtools, so a test can compare cmuts hmm against an independent count.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass

from programs import samtools, samtools_into

# The columns of a SAM record, counting from zero.
FLAG_COLUMN = 1
REFERENCE_COLUMN = 2
MAPQ_COLUMN = 4
CIGAR_COLUMN = 5
SEQUENCE_COLUMN = 9
FIRST_TAG_COLUMN = 11

# The flags samtools is given to select a strand.
STRAND_FLAGS = {
    "forward,reverse": (),
    "forward": ("-F", "16"),
    "reverse": ("-f", "16"),
}

# The SAM spec's "unavailable" mapping quality. samtools accepts a read with
# it at every -q, and cmuts hmm rejects one at every threshold.
UNAVAILABLE_MAPQ = 255

UNMAPPED_FLAG = ("-f", "4")
MAPPED_FLAG = ("-F", "4")

# Unmapped and secondary reads, which cmuts hmm rejects under every criterion.
NOT_COUNTED_FLAGS = ("-F", "0x104")


@dataclass(frozen=True)
class Record:
    """One alignment, holding the columns a test reads."""

    flag: int
    reference: str
    mapq: int
    cigar: str
    sequence: str
    tags: tuple


def _parse(line: str) -> Record:
    columns = line.split("\t")

    return Record(
        flag=int(columns[FLAG_COLUMN]),
        reference=columns[REFERENCE_COLUMN],
        mapq=int(columns[MAPQ_COLUMN]),
        cigar=columns[CIGAR_COLUMN],
        sequence=columns[SEQUENCE_COLUMN],
        tags=tuple(columns[FIRST_TAG_COLUMN:]),
    )


def record_lines(bam, *flags) -> list:
    """Returns every alignment of a file as the line samtools prints for it."""
    return [line for line in samtools("view", *flags, bam).splitlines() if line]


def records(bam, *flags) -> list:
    """Returns every alignment of a file, parsed."""
    return [_parse(line) for line in record_lines(bam, *flags)]


def count_records(bam, *flags) -> int:
    """Returns how many alignments pass a set of samtools flags."""
    return int(samtools("view", "-c", *flags, bam))


def header_text(bam) -> str:
    """Returns the header samtools prints for a file."""
    return samtools("view", "-H", bam)


def header_lines(bam, prefix: str) -> list:
    """Returns the header lines starting with a tag."""
    return [line for line in header_text(bam).splitlines() if line.startswith(prefix)]


def references_with_reads(bam) -> set:
    """Returns the names of the references that a mapped read aligned to."""
    return {record.reference for record in records(bam, *MAPPED_FLAG)}


# ---------------------------------------------------------------------------
# Reading a FASTA
# ---------------------------------------------------------------------------


def sequences(fasta) -> dict:
    """Returns every record of a FASTA, keyed by name, in file order."""
    records, name = {}, None

    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name:
            records[name].append(line.strip())

    return {name: "".join(parts) for name, parts in records.items()}


def reference_lengths(fasta) -> dict:
    """Returns the length of every reference, keyed by name."""
    return {name: len(seq) for name, seq in sequences(fasta).items()}


def rows_by_name(fasta) -> dict:
    """Returns the row of the output each reference is written to.

    Rows are identified by position and not by name: the FASTA and the header
    declare the references in the same order, so a reference's row is where it
    sits in the FASTA.
    """
    return {name: i for i, name in enumerate(sequences(fasta))}


# ---------------------------------------------------------------------------
# The reads that pass a set of criteria
# ---------------------------------------------------------------------------


def surviving_records(data, min_mapq: int, strand: str) -> list:
    """Returns the records passing a mapping quality and a strand.

    Unmapped and secondary reads are excluded, as are the reads of unavailable
    mapping quality. cmuts hmm rejects all three under every criterion, so the
    divergence from samtools is held here and not in each test.
    """
    flags = (*NOT_COUNTED_FLAGS, "-q", str(min_mapq), *STRAND_FLAGS[strand])

    return [record for record in records(data.bam, *flags)
            if record.mapq != UNAVAILABLE_MAPQ]


def _within_length(stored: int, min_length: int, max_length: int) -> bool:
    """Returns whether a stored length falls inside the bounds. A bound of zero
    is left unapplied, matching an unset bound in cmuts hmm."""
    return ((min_length == 0 or stored >= min_length)
            and (max_length == 0 or stored <= max_length))


def samtools_kept(
    data, min_mapq: int = 0, strand: str = "forward,reverse",
    min_length: int = 0, max_length: int = 0,
) -> int:
    """Returns how many reads pass a set of criteria.

    samtools takes a strand and a mapping quality as flags but not a length, so
    the length is measured from the sequence column here.
    """
    return sum(
        1 for record in surviving_records(data, min_mapq, strand)
        if _within_length(len(record.sequence), min_length, max_length)
    )


def samtools_length_histogram(
    data, min_mapq: int = 0, strand: str = "forward,reverse",
) -> dict:
    """Returns the stored length of every passing read, counted per reference.

    `samtools stats` reports one histogram for the file as a whole, and its -q
    is a trimming parameter and not a bound on mapping quality, so the counts
    are taken from the records here.
    """
    counts = defaultdict(Counter)

    for record in surviving_records(data, min_mapq, strand):
        counts[record.reference][len(record.sequence)] += 1

    return counts


def assert_counts_agree(summary, data, criteria: dict):
    """Asserts every claim a run makes about how many reads it saw, against
    samtools and against the run's own totals."""
    assert summary.kept == samtools_kept(data, **criteria), "surviving reads"

    # Every mapped read is either kept or rejected.
    assert summary.kept + summary.rejected == data.mapped, "reads accounted for"

    # An unmapped read aligns nowhere, so no filter can reach it.
    assert summary.unmapped == data.unmapped, "unmapped reads"

    # A reference that received any mapped read gets a row, whether or not any
    # of those reads passed the filter.
    assert summary.rows == data.touched, "references written"


# ---------------------------------------------------------------------------
# The tags samtools computes
# ---------------------------------------------------------------------------


def md_and_nm_tags(bam) -> list:
    """Returns the MD and NM tags of every record, in file order."""
    return [
        tuple(tag for tag in record.tags if tag.startswith(("MD:", "NM:")))
        for record in records(bam)
    ]


def recomputed_md_and_nm_tags(data, directory) -> list:
    """Returns the MD and NM tags samtools computes from the alignment and the
    reference alone."""
    output = directory / "calmd.bam"

    samtools("faidx", data.fasta)
    samtools_into(output, "calmd", "-b", data.bam, data.fasta)

    return md_and_nm_tags(output)


def reference_span(cigar: str) -> int:
    """Returns how many reference positions a CIGAR consumes."""
    span, digits = 0, ""

    for character in cigar:
        if character.isdigit():
            digits += character
        else:
            if character in "MDN=X":
                span += int(digits)
            digits = ""

    return span
