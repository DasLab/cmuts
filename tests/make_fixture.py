"""Generate the synthetic alignments the filtering tests run against.

Writes a SAM and a matching FASTA. The reads are a deterministic cross of
mapping quality, strand and length, so that every filter has cases on both
sides of it, plus the two reference shapes that are easy to get wrong: one
whose reads are all rejected by a modest threshold, and one with no reads at
all.

Nothing here depends on what the processing step computes. The tests compare
how many reads survived, which samtools can be asked independently.
"""

import argparse
import pathlib

# Reference name, length. Lengths differ so rows are padded unevenly.
REFERENCES = [
    ("ref_short", 60),
    ("ref_mid", 150),
    ("ref_long", 400),
    ("ref_lowqual", 100),   # every read below any useful MAPQ threshold
    ("ref_empty", 200),     # no reads at all
]

MAPQS = [0, 1, 5, 30, 60, 255]
STRANDS = [0, 16]           # forward, reverse
LENGTHS = [20, 35, 50, 90, 150]

BASES = "ACGT"


def sequence_for(name, length):
    """A stable, arbitrary sequence; content is irrelevant to read counting."""
    return "".join(BASES[(i * 7 + len(name)) % 4] for i in range(length))


def reads_for(reference, reflen, sequence, start_id):
    """One read per (mapq, strand, length) that fits inside the reference."""
    lengths = [n for n in LENGTHS if n <= reflen]
    mapqs = [0] if reference == "ref_lowqual" else MAPQS
    rows, n = [], start_id

    for mapq in mapqs:
        for flag in STRANDS:
            for length in lengths:
                pos = 1 + (n * 3) % (reflen - length + 1)
                rows.append((
                    f"read{n}", flag, reference, pos, mapq, f"{length}M",
                    sequence[pos - 1:pos - 1 + length], "I" * length,
                ))
                n += 1

    return rows, n


def write_fixture(directory):
    directory.mkdir(parents=True, exist_ok=True)
    sequences = {name: sequence_for(name, length) for name, length in REFERENCES}

    with open(directory / "fixture.fasta", "w") as fasta:
        for name, length in REFERENCES:
            fasta.write(f">{name}\n{sequences[name]}\n")

    rows, next_id = [], 0
    for name, length in REFERENCES:
        if name == "ref_empty":
            continue
        produced, next_id = reads_for(name, length, sequences[name], next_id)
        rows.extend(produced)

    with open(directory / "fixture.sam", "w") as sam:
        sam.write("@HD\tVN:1.6\tSO:unsorted\n")
        for name, length in REFERENCES:
            sam.write(f"@SQ\tSN:{name}\tLN:{length}\n")

        for name, flag, ref, pos, mapq, cigar, seq, qual in rows:
            sam.write(f"{name}\t{flag}\t{ref}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t"
                      f"{seq}\t{qual}\tNM:i:0\tMD:Z:{len(seq)}\n")

        # Unmapped reads belong to no reference and must be counted apart.
        for i in range(7):
            sam.write(f"unmapped{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n")

    return len(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    args = parser.parse_args()

    mapped = write_fixture(args.directory)
    print(f"{len(REFERENCES)} references, {mapped} mapped reads, 7 unmapped")
