"""Running the programs under test, and asking samtools what to expect.

samtools is the oracle: it decides independently how many reads should survive
a set of criteria, and the tests ask only whether cmuts agrees. Nothing here
inspects what the processing step computed, so these tests stay valid however
that changes.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass, replace
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CMUTS = ROOT / "build" / "cmuts"
CMUTS_GEN = ROOT / "build" / "cmuts-gen"

STRAND_FLAGS = {"both": (), "forward": ("-F", "16"), "reverse": ("-f", "16")}


def _run(command, **kwargs):
    return subprocess.run(
        [str(word) for word in command], check=True, capture_output=True, text=True, **kwargs
    )


def _option(name: str) -> str:
    return "--" + name.replace("_", "-")


# ---------------------------------------------------------------------------
# Generated data
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Dataset:
    """Alignments and the reference they came from, with the totals that hold
    whatever filter is applied afterwards.

    The alignments may be spread over several files, which cmuts reads as one.
    The totals are of them all, so they carry over a split.
    """

    bams: tuple
    fasta: Path
    mapped: int
    unmapped: int
    touched: int

    @property
    def bam(self) -> Path:
        """The one file, for the oracles that read a single one. Raises where
        the alignments are spread over several."""
        (bam,) = self.bams
        return bam


def generate(directory, name: str, **parameters) -> Dataset:
    prefix = Path(directory) / name
    command = [CMUTS_GEN, "-o", prefix]

    for key, value in parameters.items():
        command += [_option(key), value]

    _run(command)

    bam = Path(f"{prefix}.bam")
    return Dataset(
        bams=(bam,),
        fasta=Path(f"{prefix}.fasta"),
        mapped=_count(bam, "-F", 4),
        unmapped=_count(bam, "-f", 4),
        touched=len(references_with_reads(bam)),
    )


def converted(data: Dataset, directory, fmt: str) -> Dataset:
    """The same alignments in another format.

    CRAM stores sequence as differences from a reference, so it needs one to be
    written as well as read; the totals are unchanged by the conversion.
    """
    if fmt == "bam":
        return data

    output = Path(directory) / f"converted.{fmt}"
    flags = ("-T", str(data.fasta)) if fmt == "cram" else ()

    _run(["samtools", "view", "-h", "-O", fmt.upper(), *flags, "-o", output, data.bam])

    return replace(data, bams=(output,))


def dealt_out(data: Dataset, directory, parts: int) -> Dataset:
    """The same alignments dealt between several files.

    Records keep their relative order, so every part is coordinate sorted and
    holds an interleaved share of the references rather than a run of them.
    Nothing is added or lost, so the totals carry over.
    """
    header = _run(["samtools", "view", "-H", data.bam]).stdout
    records = _records(data.bam)
    written = []

    for i in range(parts):
        sam = Path(directory) / f"part{i}.sam"
        sam.write_text(header + "".join(line + "\n" for line in records[i::parts]))

        bam = Path(directory) / f"part{i}.bam"
        _run(["samtools", "view", "-b", "-o", bam, sam])
        written.append(bam)

    return replace(data, bams=tuple(written))


def reheadered(data: Dataset, directory, transform) -> Dataset:
    """The same alignments behind a header the transform has rewritten.

    Only the header changes, so the totals carry over and any difference in
    what cmuts does is down to what the header says.
    """
    header = Path(directory) / "header.sam"
    header.write_text(transform(_run(["samtools", "view", "-H", data.bam]).stdout))

    bam = Path(directory) / "reheadered.bam"
    with open(bam, "wb") as handle:
        subprocess.run(["samtools", "reheader", str(header), str(data.bam)],
                       check=True, stdout=handle, stderr=subprocess.DEVNULL)

    return replace(data, bams=(bam,))


# ---------------------------------------------------------------------------
# What samtools says
# ---------------------------------------------------------------------------


def _count(bam, *flags) -> int:
    return int(_run(["samtools", "view", "-c", *flags, bam]).stdout)


def _records(bam, *flags):
    return [line for line in _run(["samtools", "view", *flags, bam]).stdout.splitlines() if line]


def references_with_reads(bam) -> set:
    return {line.split("\t")[2] for line in _records(bam, "-F", "4")}


def samtools_kept(
    data: Dataset, min_mapq: int = 0, strand: str = "both",
    min_length: int = 0, max_length: int = 0,
) -> int:
    """Reads surviving a set of criteria.

    Strand and mapping quality are flags samtools knows; length it does not, so
    the sequence column is measured directly. A bound of zero is not applied,
    which is what cmuts does with one left unset.
    """
    flags = ("-F", "4", "-q", str(min_mapq), *STRAND_FLAGS[strand])
    lengths = (len(line.split("\t")[9]) for line in _records(data.bam, *flags))

    return sum(
        1
        for n in lengths
        if (min_length == 0 or n >= min_length) and (max_length == 0 or n <= max_length)
    )


def md_and_nm_tags(bam):
    """The MD and NM of every record, in file order."""
    return [
        tuple(tag for tag in line.split("\t")[11:] if tag.startswith(("MD:", "NM:")))
        for line in _records(bam)
    ]


def recomputed_md_and_nm_tags(data: Dataset, directory):
    """What samtools makes of MD and NM given only the alignment and the
    reference, which is independent of how the generator wrote them."""
    output = Path(directory) / "calmd.bam"

    _run(["samtools", "faidx", data.fasta])
    with open(output, "wb") as handle:
        subprocess.run(
            ["samtools", "calmd", "-b", str(data.bam), str(data.fasta)],
            check=True, stdout=handle, stderr=subprocess.DEVNULL,
        )

    return md_and_nm_tags(output)


# ---------------------------------------------------------------------------
# What cmuts says
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Summary:
    kept: int
    rejected: int
    rows: int
    unmapped: int


def read_summary(path) -> Summary:
    with h5py.File(path, "r") as output:
        reads = output["reads"][:]
        rejected = output["reads_filtered"][:]

        return Summary(
            kept=int(np.nansum(reads, dtype=np.float64)),
            rejected=int(np.nansum(rejected, dtype=np.float64)),
            rows=int((~np.isnan(reads)).sum()),
            unmapped=int(output.attrs["reads_unmapped"]),
        )


def _cmuts_command(data: Dataset, output, workers: int, options):
    command = [CMUTS, "-f", data.fasta, "-o", output, "-j", workers]

    for key, value in options.items():
        # A flag carries no value of its own.
        command += [_option(key)] if value is True else [_option(key), value]

    return [*command, *data.bams]


def run_cmuts(data: Dataset, output, workers: int = 4, **options) -> Summary:
    _run(_cmuts_command(data, output, workers, options))
    return read_summary(output)


def try_cmuts(data: Dataset, output, workers: int = 4, **options):
    """Runs it whether or not it succeeds, for the paths that should fail."""
    return subprocess.run(
        [str(word) for word in _cmuts_command(data, output, workers, options)],
        capture_output=True, text=True,
    )


def _datasets_agree(a, b) -> bool:
    # NaN marks a reference no read reached, so two outputs agree where both
    # hold one. Only the counting datasets can carry it; the names cannot.
    return np.array_equal(a[:], b[:], equal_nan=np.issubdtype(a.dtype, np.floating))


def outputs_agree(first, second) -> bool:
    """Whether two outputs hold the same thing.

    Over every dataset they carry, rather than a list of them kept here: a
    quantity added to the accumulator would otherwise be written by the code and
    checked by nothing.
    """
    with h5py.File(first, "r") as a, h5py.File(second, "r") as b:
        return set(a) == set(b) and all(_datasets_agree(a[k], b[k]) for k in a)


def counted_fields(path):
    """The datasets holding counts, which are the ones arithmetic applies to."""
    with h5py.File(path, "r") as output:
        return [k for k in output if np.issubdtype(output[k].dtype, np.floating)]
