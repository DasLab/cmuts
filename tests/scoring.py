"""Tools for building what cmuts score reads, and for computing what it should
write.

A row of the table depends on the reactivities, the sequence and the pairing
alone, so the files are written by hand. Each metric is derived here from its
definition, over every pair of positions, so that the oracle and the program do
not arrive at an answer the same way.
"""

from __future__ import annotations

import numpy as np

from outputs import COVERAGE, REACTIVITY, field_of

UNPAIRED = "."
BRACKETS = "()[]{}<>"

# The columns of the table, in the order they are written.
COLUMNS = ("reference", "paired", "unpaired", "auroc", "auprc",
           "mean_paired", "mean_unpaired")

# The bases a run scores when --bases is not given.
ALL_BASES = "ACGU"


# ---------------------------------------------------------------------------
# Writing what a run reads
# ---------------------------------------------------------------------------


def write_fasta(path, sequences: dict):
    """Writes the references in the order given, which is the order the rows of
    the file follow."""
    with open(path, "w") as out:
        for name, sequence in sequences.items():
            out.write(f">{name}\n{sequence}\n")

    return path


def write_structures(path, pairings: dict, sequences: dict = None):
    """Writes dot bracket records. A record carries its sequence before the
    pairing where one is given, which a reader must pass over."""
    with open(path, "w") as out:
        for name, pairing in pairings.items():
            out.write(f">{name}\n")

            if sequences:
                out.write(f"{sequences[name]}\n")

            out.write(f"{pairing}\n")

    return path


def alternating(length: int) -> str:
    """A pairing leaving every third base open, so that both classes are held
    whatever the length."""
    return "".join(UNPAIRED if i % 3 == 0 else "(" for i in range(length))


def sequences_of(path) -> dict:
    """Reads a FASTA by record name, for a structure built to match one that a
    run generated."""
    sequences, name = {}, None

    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = ""
        elif name is not None:
            sequences[name] += line.strip()

    return sequences


# ---------------------------------------------------------------------------
# Reading the table
# ---------------------------------------------------------------------------


def rows_of(text: str) -> dict:
    """Returns the table by reference name, each row keyed by column. The
    header is checked here, so a test reads a value by name and not by
    position."""
    lines = text.splitlines()

    assert lines, "the table is empty"
    assert lines[0] == ",".join(COLUMNS)

    rows = {}

    for line in lines[1:]:
        cells = line.split(",")

        assert len(cells) == len(COLUMNS), line

        rows[cells[0]] = {name: float(cell)
                          for name, cell in zip(COLUMNS[1:], cells[1:])}

    return rows


def order_of(text: str) -> list:
    """Returns the references the table names, in the order it names them."""
    return [line.split(",")[0] for line in text.splitlines()[1:]]


# ---------------------------------------------------------------------------
# The oracle
# ---------------------------------------------------------------------------


def kept(path, tid, sequence, pairing, bases=ALL_BASES, min_coverage=0.0):
    """Returns the value of every position a run scores, and whether the
    structure leaves each open.

    A position is kept where the structure names it, the base is one of those
    asked for, the reactivity is finite, and the coverage clears the floor.
    """
    rates = np.asarray(field_of(path, REACTIVITY))[tid]
    coverage = np.asarray(field_of(path, COVERAGE))[tid]

    values, opened = [], []

    for i, (base, mark) in enumerate(zip(sequence, pairing)):
        named = mark == UNPAIRED or mark in BRACKETS
        wanted = base.upper().replace("T", "U") in bases

        if not named or not wanted:
            continue

        if not np.isfinite(rates[i]) or coverage[i] < min_coverage:
            continue

        values.append(float(rates[i]))
        opened.append(mark == UNPAIRED)

    return np.array(values), np.array(opened, dtype=bool)


def auroc(values, opened) -> float:
    """The chance an unpaired base outranks a paired one, counted over every
    pair of them, with a tie counted as half."""
    unpaired, paired = values[opened], values[~opened]

    wins = int((unpaired[:, None] > paired[None, :]).sum())
    ties = int((unpaired[:, None] == paired[None, :]).sum())

    return (wins + (ties / 2)) / (unpaired.size * paired.size)


def auprc(values, opened) -> float:
    """The average precision, with the unpaired bases as the positive class and
    each distinct value in turn as the threshold."""
    total = int(opened.sum())
    previous, average = 0.0, 0.0

    for threshold in sorted(set(values.tolist()), reverse=True):
        at_or_above = values >= threshold
        found = int((at_or_above & opened).sum())
        recall = found / total

        average += (recall - previous) * (found / int(at_or_above.sum()))
        previous = recall

    return average


def expected(values, opened) -> dict:
    """Everything one row of the table should hold."""
    return {
        "paired": float((~opened).sum()),
        "unpaired": float(opened.sum()),
        "auroc": auroc(values, opened),
        "auprc": auprc(values, opened),
        "mean_paired": float(values[~opened].mean()),
        "mean_unpaired": float(values[opened].mean()),
    }
