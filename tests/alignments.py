"""Tools for generating alignments and transforming them.

A transform returns a new dataset that differs from its input in one respect,
so that a test can attribute any change in the result to that one difference.
The totals carry over by construction, and are not counted a second time.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

from oracle import (
    FLAG_COLUMN,
    MAPPED_FLAG,
    UNMAPPED_FLAG,
    count_records,
    header_text,
    record_lines,
    references_with_reads,
    sequences,
)
from programs import run_generator, samtools, samtools_into

# The bits of the SAM flag a transform sets or reads.
PAIRED_BIT = 0x1
SECONDARY_BIT = 0x100
UNMAPPED_BIT = 0x4

BAM = "bam"
SAM = "sam"
CRAM = "cram"

# The format cmuts gen writes, which the others are converted from.
NATIVE = BAM

FORMATS = (BAM, SAM, CRAM)

# The formats that store a read's sequence in the file. A CRAM stores it as
# differences from a reference, so it cannot be decoded against any other
# sequence and a test giving it one measures htslib.
SELF_CONTAINED = (BAM, SAM)


@dataclass(frozen=True)
class Dataset:
    """Alignments, the reference they were made against, and the totals that no
    filter can change.

    The alignments may be spread over several files, which cmuts hmm reads as
    one, so the totals are of them all. Every file is in the format fmt names,
    which a transform preserves.
    """

    name: str
    bams: tuple
    fasta: Path
    mapped: int
    unmapped: int
    touched: int
    fmt: str = NATIVE

    @property
    def bam(self) -> Path:
        """The single alignment file. Raises where the alignments are spread
        over several."""
        (bam,) = self.bams
        return bam


def measure_dataset(name: str, bams, fasta) -> Dataset:
    """Builds a dataset whose totals are measured from the files, so that a
    caller cannot declare a total the files do not support."""
    bams = tuple(bams)

    return Dataset(
        name=name,
        bams=bams,
        fasta=fasta,
        mapped=sum(count_records(bam, *MAPPED_FLAG) for bam in bams),
        unmapped=sum(count_records(bam, *UNMAPPED_FLAG) for bam in bams),
        touched=len({reference for bam in bams
                     for reference in references_with_reads(bam)}),
    )


def generate(directory, name: str, **parameters) -> Dataset:
    """Generates a dataset with cmuts gen."""
    prefix = Path(directory) / name
    run_generator(prefix, parameters)

    return measure_dataset(name, (Path(f"{prefix}.bam"),), Path(f"{prefix}.fasta"))


# ---------------------------------------------------------------------------
# Rewriting the alignments
# ---------------------------------------------------------------------------


def _rewrite(source: Path, output: Path, fmt: str, fasta: Path) -> Path:
    """Writes an alignment file out in the named format.

    CRAM stores a sequence as its differences from a reference, so writing one
    requires a reference as well as reading one does.
    """
    flags = ("-T", str(fasta)) if fmt == CRAM else ()

    samtools("view", "-h", "-O", fmt.upper(), *flags, "-o", output, source)

    return output


def _write_sam(path: Path, header: str, lines) -> Path:
    path.write_text(header + "".join(line + "\n" for line in lines))
    return path


def _rebuild(data: Dataset, directory, stem: str, header: str, lines) -> Path:
    """Writes records out in the format the dataset is stored in, so that a
    transform leaves the format it was given.

    Records are written as SAM first, that being the format they are held in
    here. A dataset stored as SAM is finished at that point.
    """
    sam = _write_sam(Path(directory) / f"{stem}.sam", header, lines)

    if data.fmt == SAM:
        return sam

    return _rewrite(sam, Path(directory) / f"{stem}.{data.fmt}", data.fmt, data.fasta)


def convert_format(data: Dataset, directory, fmt: str) -> Dataset:
    """Converts the alignments to another format."""
    if fmt == data.fmt:
        return data

    output = _rewrite(data.bam, Path(directory) / f"{data.name}.{fmt}", fmt, data.fasta)

    return replace(data, bams=(output,), fmt=fmt)


def split_across_files(data: Dataset, directory, parts: int) -> Dataset:
    """Deals the alignments between several files.

    Records keep their relative order, so every part is coordinate sorted and
    holds an interleaved share of the references.
    """
    header = header_text(data.bam)
    lines = record_lines(data.bam)

    written = [_rebuild(data, directory, f"part{i}", header, lines[i::parts])
               for i in range(parts)]

    return replace(data, bams=tuple(written))


def _with_flag_bit(line: str, bit: int) -> str:
    columns = line.split("\t")
    columns[FLAG_COLUMN] = str(int(columns[FLAG_COLUMN]) | bit)

    return "\t".join(columns)


def _is_mapped(line: str) -> bool:
    return not int(line.split("\t")[FLAG_COLUMN]) & UNMAPPED_BIT


def mark_secondary(data: Dataset, directory, every: int):
    """Marks a share of the mapped reads as secondary alignments.

    Only the flag changes, so the totals and the sort order carry over. Returns
    the dataset and how many reads were marked.
    """
    lines, marked = [], 0

    for i, line in enumerate(record_lines(data.bam)):
        if _is_mapped(line) and i % every == 0:
            line = _with_flag_bit(line, SECONDARY_BIT)
            marked += 1

        lines.append(line)

    rebuilt = _rebuild(data, directory, "secondary", header_text(data.bam), lines)

    return replace(data, bams=(rebuilt,)), marked


def mark_paired(data: Dataset, directory) -> Dataset:
    """Marks every read as one mate of a pair.

    Only the flag changes, so the totals and the sort order carry over.
    """
    lines = [_with_flag_bit(line, PAIRED_BIT) for line in record_lines(data.bam)]

    rebuilt = _rebuild(data, directory, "paired", header_text(data.bam), lines)

    return replace(data, bams=(rebuilt,))


# ---------------------------------------------------------------------------
# Rewriting the reference
# ---------------------------------------------------------------------------

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def write_fasta(records: dict, path: Path) -> Path:
    """Writes a FASTA holding the given sequences, by name, in the order
    given."""
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()))

    return path


def replace_bases(data: Dataset, directory, only=None) -> Dataset:
    """Complements the bases of the reference, leaving every name and length
    unchanged.

    A name and a length therefore describe the new reference as well as they
    described the old one, and only the bases distinguish the two.
    """
    records = {
        name: seq.translate(COMPLEMENT) if only is None or name in only else seq
        for name, seq in sequences(data.fasta).items()
    }

    return replace(data, fasta=write_fasta(records, Path(directory) / "substituted.fasta"))


def lengthen_references(data: Dataset, directory, extra: int = 500) -> Dataset:
    """Extends every reference of the FASTA. The headers still declare the
    original lengths, so every record is longer than its declaration."""
    records = {
        name: seq + ("A" * extra) for name, seq in sequences(data.fasta).items()
    }

    return replace(data, fasta=write_fasta(records, Path(directory) / "lengthened.fasta"))


def rename_references(data: Dataset, directory) -> Dataset:
    """Renames every reference of the FASTA, leaving the bases and their order
    unchanged."""
    records = {f"contig{name}": seq for name, seq in sequences(data.fasta).items()}

    return replace(data, fasta=write_fasta(records, Path(directory) / "renamed.fasta"))


# ---------------------------------------------------------------------------
# Rewriting the header
# ---------------------------------------------------------------------------


def replace_header(data: Dataset, directory, transform) -> Dataset:
    """Rewrites the header of the alignments with a transform over its text.

    Only the header changes, so the totals carry over and any difference in
    what cmuts hmm does follows from the header contents.

    samtools reheader reads a BAM or a CRAM. A SAM is rewritten here instead,
    which leaves the records untouched in the same way.
    """
    header = transform(header_text(data.bam))

    if data.fmt == SAM:
        rewritten = _write_sam(Path(directory) / "reheadered.sam",
                               header, record_lines(data.bam))
    else:
        source = Path(directory) / "header.sam"
        source.write_text(header)
        rewritten = samtools_into(Path(directory) / f"reheadered.{data.fmt}",
                                  "reheader", source, data.bam)

    return replace(data, bams=(rewritten,))


def _field_of(fields: list, tag: str):
    return next((field[len(tag):] for field in fields if field.startswith(tag)), None)


def replace_checksums(data: Dataset, directory, checksum) -> Dataset:
    """Replaces the M5 of each @SQ line with the value checksum returns, and
    removes the M5 where checksum returns None.

    checksum is passed the reference name and the M5 that cmuts gen wrote.
    cmuts gen writes a correct M5 for every reference, so a test alters what is
    already in the header and computes no digest of its own.
    """
    def replace_in_line(line):
        fields = line.split("\t")
        kept = [field for field in fields if not field.startswith("M5:")]
        m5 = checksum(_field_of(fields, "SN:"), _field_of(fields, "M5:"))

        return "\t".join(kept + ([f"M5:{m5}"] if m5 is not None else []))

    def transform(text):
        return "".join(
            (replace_in_line(line) if line.startswith("@SQ") else line) + "\n"
            for line in text.splitlines()
        )

    return replace_header(data, directory, transform)


# ---------------------------------------------------------------------------
# Alignments written by hand
# ---------------------------------------------------------------------------

# Clear of every bound the filter applies by default, so a hand-built read is
# dropped only where a test sets a bound that excludes it.
MAPQ = 60

# One score for every base, high enough that a misread does not explain a
# difference from the reference.
BASE_QUALITY = "I"


def build_placements(directory, name, reference, read, cigars) -> Dataset:
    """Writes the same read against the same reference under each CIGAR, one
    reference per CIGAR, so that a single run scores them all.

    Every reference holds the same sequence and one read, so the rows
    of the output differ only in how the alignment was written. Each reference
    takes the name of the CIGAR it holds, which a failure reports.

    The file is written as SAM, which cmuts hmm reads as it reads BAM. The
    totals are declared here because this writes every record itself.
    """
    prefix = Path(directory) / name
    fasta = Path(f"{prefix}.fasta")

    fasta.write_text("".join(f">{cigar}\n{reference}\n" for cigar in cigars))

    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    lines += [f"@SQ\tSN:{cigar}\tLN:{len(reference)}" for cigar in cigars]
    lines += [
        f"r\t0\t{cigar}\t1\t{MAPQ}\t{cigar}\t*\t0\t0\t{read}"
        f"\t{BASE_QUALITY * len(read)}"
        for cigar in cigars
    ]

    sam = _write_sam(Path(f"{prefix}.sam"), "", lines)

    return Dataset(
        name=name, bams=(sam,), fasta=fasta,
        mapped=len(cigars), unmapped=0, touched=len(cigars),
    )
