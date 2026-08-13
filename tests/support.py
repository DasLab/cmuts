"""Running the programs under test, and asking samtools what to expect.

samtools is the oracle: it decides independently how many reads should survive
a set of criteria, and the tests ask only whether cmuts agrees. Nothing here
inspects what the processing step computed, so these tests stay valid however
that changes.

Alignments built here by hand have no oracle. They carry an ambiguity the
reference cannot resolve, and what is asked of them is that cmuts answers it the
same way however the ambiguity was written down.
"""

from __future__ import annotations

import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path

import h5py
import numpy as np

from outputs import RATES

ROOT = Path(__file__).resolve().parent.parent

CMUTS = "cmuts"
CMUTS_GEN = "cmuts-gen"
CMUTS_SUB = "cmuts-sub"
PROGRAMS = (CMUTS, CMUTS_GEN, CMUTS_SUB)


def located(name: str) -> Path | None:
    """Finds a program on PATH, disregarding any copy from outside this
    repository.

    make check puts the build directory first on PATH, which is what decides
    between an ordinary build and a sanitized one. A search that got as far as
    an installed cmuts would test whatever was installed, whenever that was, so
    anything outside the repository counts as not found.
    """
    found = shutil.which(name)

    if found is None:
        return None

    path = Path(found).resolve()

    return path if ROOT in path.parents else None


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
    """Alignments and the reference they came from, with the totals that no
    later filter can change.

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


def counted(bams, fasta) -> Dataset:
    """A dataset whose totals are measured from the files, so that a caller
    cannot declare a total the files do not support.

    The transforms below carry their totals over by construction and keep them
    instead of counting again, which on a CRAM is not cheap.
    """
    bams = tuple(bams)

    return Dataset(
        bams=bams,
        fasta=fasta,
        mapped=sum(_count(bam, "-F", 4) for bam in bams),
        unmapped=sum(_count(bam, "-f", 4) for bam in bams),
        touched=len({name for bam in bams for name in references_with_reads(bam)}),
    )


def generate(directory, name: str, **parameters) -> Dataset:
    prefix = Path(directory) / name
    command = [CMUTS_GEN, "-o", prefix]

    for key, value in parameters.items():
        command += [_option(key), value]

    _run(command)

    return counted((Path(f"{prefix}.bam"),), Path(f"{prefix}.fasta"))


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
    lines = records(data.bam)
    written = []

    for i in range(parts):
        sam = Path(directory) / f"part{i}.sam"
        sam.write_text(header + "".join(line + "\n" for line in lines[i::parts]))

        bam = Path(directory) / f"part{i}.bam"
        _run(["samtools", "view", "-b", "-o", bam, sam])
        written.append(bam)

    return replace(data, bams=tuple(written))


def with_secondary(data: Dataset, directory, every: int):
    """The same alignments with a share of the mapped ones marked secondary.

    Only the flag changes, so the totals and the sort order carry over. Returns
    the count marked, that being what the filter is expected to remove.
    """
    header = _run(["samtools", "view", "-H", data.bam]).stdout
    marked = 0
    lines = []

    for i, record in enumerate(records(data.bam)):
        columns = record.split("\t")
        flag = int(columns[1])

        if not flag & 0x4 and i % every == 0:
            columns[1] = str(flag | 0x100)
            marked += 1

        lines.append("\t".join(columns))

    sam = Path(directory) / "secondary.sam"
    sam.write_text(header + "".join(line + "\n" for line in lines))

    bam = Path(directory) / "secondary.bam"
    _run(["samtools", "view", "-b", "-o", bam, sam])

    return replace(data, bams=(bam,)), marked


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def written(records, path):
    """A FASTA holding the given sequences, by name, in the order given."""
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()))

    return path


def substituted(data: Dataset, directory, only=None) -> Dataset:
    """The same alignments against a FASTA agreeing on every name and length
    and holding different bases.

    Complemented, so a name and a length describe the reference as well as they
    did and nothing but the bases can tell the two apart.
    """
    records = {
        name: seq.translate(COMPLEMENT) if only is None or name in only else seq
        for name, seq in sequences(data.fasta).items()
    }

    return replace(data, fasta=written(records, Path(directory) / "substituted.fasta"))


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
# Hand-built alignments
# ---------------------------------------------------------------------------

# Clear of every bound the filter applies by default, so a hand-built read is
# dropped only where a test asks for it.
MAPQ = 60

# One score for every base, high enough that a misread is not what explains a
# difference from the reference.
BASE_QUALITY = "I"


def placements(directory, name, reference, read, cigars) -> Dataset:
    """The same read against the same reference under each CIGAR, one reference
    per CIGAR so that a single run scores them all.

    Every reference holds the same sequence and carries one read, so the rows of
    the output differ only in how the alignment was written. Each takes the name of
    the CIGAR it carries, which is what a failure has to report.

    SAM is written directly: cmuts reads it as it reads BAM, and neither an
    index nor a conversion is needed for a file this size.

    The totals are declared rather than counted, as the ones above are: this
    writes every record itself, so what the file holds is not in question.
    """
    prefix = Path(directory) / name
    fasta = Path(f"{prefix}.fasta")
    sam = Path(f"{prefix}.sam")

    fasta.write_text("".join(f">{cigar}\n{reference}\n" for cigar in cigars))

    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    lines += [f"@SQ\tSN:{cigar}\tLN:{len(reference)}" for cigar in cigars]
    lines += [
        f"r\t0\t{cigar}\t1\t{MAPQ}\t{cigar}\t*\t0\t0\t{read}"
        f"\t{BASE_QUALITY * len(read)}"
        for cigar in cigars
    ]

    sam.write_text("".join(line + "\n" for line in lines))

    return Dataset(
        bams=(sam,), fasta=fasta,
        mapped=len(cigars), unmapped=0, touched=len(cigars),
    )


# ---------------------------------------------------------------------------
# What samtools says
# ---------------------------------------------------------------------------


def reference_lengths(fasta) -> dict:
    """The length of every reference, by name, which is how far its own row of
    an output reaches."""
    return {name: len(seq) for name, seq in sequences(fasta).items()}


def sequences(fasta):
    """Every record of a FASTA, by name, in file order."""
    records, name = {}, None

    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name:
            records[name].append(line.strip())

    return {name: "".join(parts) for name, parts in records.items()}


def _count(bam, *flags) -> int:
    return int(_run(["samtools", "view", "-c", *flags, bam]).stdout)


def records(bam, *flags):
    """Every alignment of a file, as the lines samtools prints for it."""
    return [line for line in _run(["samtools", "view", *flags, bam]).stdout.splitlines() if line]


def references_with_reads(bam) -> set:
    return {line.split("\t")[2] for line in records(bam, "-F", "4")}


# The SAM spec's "unavailable", which is not a score above every threshold.
UNAVAILABLE_MAPQ = 255

# Where a record holds its mapping quality and its stored sequence.
MAPQ_COLUMN = 4
SEQUENCE_COLUMN = 9


def surviving(data: Dataset, min_mapq: int, strand: str) -> list:
    """The records passing a mapping quality and a strand.

    Unmapped and secondary reads are excluded, as cmuts excludes them under
    every criterion, and so are the reads of unavailable mapping quality:
    samtools admits those at every -q and cmuts refuses them at every bound,
    which is a deliberate divergence and belongs to the oracle rather than to
    each test consulting it.
    """
    flags = ("-F", "0x104", "-q", str(min_mapq), *STRAND_FLAGS[strand])

    return [
        line for line in records(data.bam, *flags)
        if int(line.split("\t")[MAPQ_COLUMN]) != UNAVAILABLE_MAPQ
    ]


def samtools_kept(
    data: Dataset, min_mapq: int = 0, strand: str = "both",
    min_length: int = 0, max_length: int = 0,
) -> int:
    """Reads surviving a set of criteria.

    samtools takes strand and mapping quality as flags but not length, so this
    measures the sequence column directly. A bound of zero is left unapplied,
    matching what cmuts does with one unset.
    """
    lengths = (
        len(line.split("\t")[SEQUENCE_COLUMN])
        for line in surviving(data, min_mapq, strand)
    )

    return sum(
        1
        for n in lengths
        if (min_length == 0 or n >= min_length) and (max_length == 0 or n <= max_length)
    )


def samtools_length_histogram(
    data: Dataset, min_mapq: int = 0, strand: str = "both",
) -> dict:
    """The stored length of every surviving read, counted per reference.

    Measured from the sequence column, as samtools_kept measures it. samtools
    has no per-reference length histogram of its own: `samtools stats` reports
    one for the file as a whole, and there its -q is a trimming parameter and
    not a bound on mapping quality, so it counts reads this filter would drop.
    """
    counts = defaultdict(Counter)

    for line in surviving(data, min_mapq, strand):
        columns = line.split("\t")
        counts[columns[2]][len(columns[SEQUENCE_COLUMN])] += 1

    return counts


def rows_by_name(fasta) -> dict:
    """Which row of an output each reference is written to.

    Rows are identified by position and not by name: the FASTA and the header must
    declare the references in the same order, so a reference's row is where it
    sits in the FASTA. Taken from the FASTA rather than from the output, that
    being the mapping a caller has to hand and so the one under test.
    """
    return {name: i for i, name in enumerate(sequences(fasta))}


def assert_counts_agree(summary, data: Dataset, criteria: dict):
    """Every claim a run makes about how many reads there were: against
    samtools, and against itself.

    Held here so that the tests of the filtering and of the formats are checking
    the same thing, rather than two lists that happen to look alike.
    """
    assert summary.kept == samtools_kept(data, **criteria), "surviving reads"

    # Every mapped read is either kept or rejected; nothing goes missing and
    # nothing is counted twice.
    assert summary.kept + summary.rejected == data.mapped, "reads accounted for"

    # Unmapped reads align nowhere, so no filter can touch them.
    assert summary.unmapped == data.unmapped, "unmapped reads"

    # A reference that received any mapped read gets a row, whether or not
    # anything survived the filter.
    assert summary.rows == data.touched, "references written"


def md_and_nm_tags(bam):
    """The MD and NM of every record, in file order."""
    return [
        tuple(tag for tag in line.split("\t")[11:] if tag.startswith(("MD:", "NM:")))
        for line in records(bam)
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
        reads = output["reads/counted"][:]
        rejected = output["reads/rejected"][:]

        return Summary(
            kept=int(np.nansum(reads, dtype=np.float64)),
            rejected=int(np.nansum(rejected, dtype=np.float64)),
            # The counts are zero-filled, so a reference the run wrote a row
            # for is one that some read reached, kept or rejected.
            rows=int(((reads + rejected) > 0).sum()),
            unmapped=int(output["reads/unmapped"][()]),
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


def _subtract_command(treated, untreated, output, options):
    command = [CMUTS_SUB, "-o", output]

    for key, value in options.items():
        # A flag carries no value of its own.
        command += [_option(key)] if value is True else [_option(key), value]

    return [*command, treated, untreated]


def run_subtract(treated, untreated, output, **options):
    """Takes the background off a run, returning where the result was written."""
    _run(_subtract_command(treated, untreated, output, options))
    return output


def try_subtract(treated, untreated, output, **options):
    """Runs it whether or not it succeeds, for the paths that should fail."""
    return subprocess.run(
        [str(word) for word in _subtract_command(treated, untreated, output, options)],
        capture_output=True, text=True,
    )


def _datasets_agree(a, b) -> bool:
    # NaN marks a reference no read reached, so two outputs agree where both
    # hold one. [()] rather than [:], a run total being a scalar that cannot be
    # sliced.
    return np.array_equal(a[()], b[()], equal_nan=np.issubdtype(a.dtype, np.floating))


def outputs_agree(first, second) -> bool:
    """Whether two outputs hold the same thing.

    Compares every dataset the files carry instead of a list kept here, so that
    adding a field to the output does not leave it unchecked.
    """
    with h5py.File(first, "r") as a, h5py.File(second, "r") as b:
        left, right = arrays_of(a), arrays_of(b)

        return set(left) == set(right) and all(
            _datasets_agree(left[k], right[k]) for k in left
        )


def arrays_of(handle) -> dict:
    """Every dataset of an open output, by path. The counts about reads sit in
    a group of their own, so a name is a path and not a key."""
    found = {}

    handle.visititems(
        lambda name, obj: found.update({name: obj})
        if isinstance(obj, h5py.Dataset) else None
    )

    return found


def counted_fields(path):
    """The datasets holding counts, which are the ones arithmetic applies to."""
    with h5py.File(path, "r") as output:
        # Counts are unsigned and coverage is a float; both are added over. The
        # rates are neither, two files holding twice the reads and one rate.
        return [k for k, d in arrays_of(output).items()
                if d.ndim >= 1 and k not in RATES]
