"""Running the programs under test, and asking samtools what to expect.

samtools is the oracle: it decides independently how many reads should survive
a set of criteria, and the tests ask only whether cmuts agrees. Nothing here
inspects what the processing step computed, so these tests stay valid however
that changes.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CMUTS = ROOT / "build" / "cmuts"
CMUTS_GEN = ROOT / "build" / "cmuts-gen"

# Datasets carrying one value per reference, which every comparison between two
# runs covers.
FIELDS = ("coverage", "mutations", "reads", "reads_filtered")

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
    whatever filter is applied afterwards."""

    bam: Path
    fasta: Path
    mapped: int
    unmapped: int
    touched: int


def generate(directory, name: str, **parameters) -> Dataset:
    prefix = Path(directory) / name
    command = [CMUTS_GEN, "-o", prefix]

    for key, value in parameters.items():
        command += [_option(key), value]

    _run(command)

    bam = Path(f"{prefix}.bam")
    return Dataset(
        bam=bam,
        fasta=Path(f"{prefix}.fasta"),
        mapped=_count(bam, "-F", 4),
        unmapped=_count(bam, "-f", 4),
        touched=len(references_with_reads(bam)),
    )


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


def run_cmuts(data: Dataset, output, workers: int = 4, **filters) -> Summary:
    command = [CMUTS, "-f", data.fasta, "-o", output, "-j", workers]

    for key, value in filters.items():
        command += [_option(key), value]

    _run([*command, data.bam])
    return read_summary(output)


def outputs_agree(first, second) -> bool:
    with h5py.File(first, "r") as a, h5py.File(second, "r") as b:
        return all(np.array_equal(a[k][:], b[k][:], equal_nan=True) for k in FIELDS)
