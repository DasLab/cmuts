"""An ambiguity the reference cannot resolve is not resolved by the aligner.

A gap inside a homopolymer can be written after any base of the run, and every
CIGAR that writes it describes the same alignment. Marginalizing over the band
must therefore reach the same arrays whichever one the aligner reported. It does
once the band is as wide as the gap, whatever the length of the run the gap sits
in: the width the reference asks for follows the gap alone.
"""

import h5py
import numpy as np
import pytest

from support import placements, run_cmuts

KINDS = ("D", "I")
GAPS = (1, 2, 3)
RUNS = (10, 30)

# Flanks pin both ends of the read. Sliding it by a base would set two of them
# against the run, which costs more than the gap that slide removes, so the
# alignment is left with nowhere to go but the run. They differ from each other
# so that neither end of the read can meet the wrong one.
LEFT = "GC"
RIGHT = "TG"

# How far two rows may stand apart and still be the same answer. A row's band is
# drawn around its own CIGAR, so the row carrying the gap is wider than the rest
# and each placement admits slightly different paths at the edges of its band.
# What that leaves is a real difference in the marginal and not rounding: it
# measures around 6e-5, three orders above the resolution the counts are stored
# at, and it shrinks as the band widens.
TOLERANCE = 1e-3

# How far apart a band too narrow for its gap must leave them. Well inside the
# smallest departure measured, which is a third of a mutation.
DIVIDED = 0.1

# Insertions carry no weight by default, which would leave every insertion case
# reading an array of zeros and agreeing with itself. What is under test is
# where an event is counted, not what a run chooses to count, so both kinds are
# given a weight here.
WEIGHTED = dict(insertion_weight=1)

CASES = [(kind, gap, run) for kind in KINDS for gap in GAPS for run in RUNS]
CONTROLS = [(kind, gap, RUNS[0]) for kind in KINDS for gap in GAPS]


def describe(case) -> str:
    kind, gap, run = case
    return f"{gap}{kind}-in-{run}"


def homopolymer(kind: str, gap: int, run: int):
    """The reference, the read, and every CIGAR placing the gap in the run.

    A deletion takes the gap out of the read and may follow any of the bases it
    leaves behind; an insertion adds it, and so has one placement more than the
    run has bases.
    """
    reference = LEFT + "A" * run + RIGHT
    read = LEFT + "A" * (run - gap if kind == "D" else run + gap) + RIGHT
    spots = range(run - gap + 1) if kind == "D" else range(run + 1)
    cigars = []

    for taken in spots:
        left = len(LEFT) + taken
        right = len(read) - left - (gap if kind == "I" else 0)
        cigars.append(f"{left}M{gap}{kind}{right}M")

    return reference, read, cigars


def written_rows(output, placements: int) -> dict:
    """Every dataset a row can be read from, by name. The reference names label
    the rows and differ by construction, so they are not among them.

    One row per placement, in every dataset. Nothing below would say otherwise:
    the rows are asserted equal to one another, so a dataset short of a row
    would be compared over fewer and agree, and the spread of one would be
    taken over fewer and be narrower.
    """
    with h5py.File(output, "r") as handle:
        rows = {name: handle[name][:] for name in handle if name != "reference"}

    for name, values in rows.items():
        assert len(values) == placements, \
            f"{name}: {len(values)} rows for {placements} placements"

    return rows


def spread(rows) -> float:
    """How far apart the rows of one dataset stand, over every position."""
    return float(np.nanmax(np.nanmax(rows, axis=0) - np.nanmin(rows, axis=0)))


def scored(kind: str, gap: int, read: str) -> int:
    """Bases of the read that meet a reference position. An inserted base
    answers to none, and so covers none."""
    return len(read) - gap if kind == "I" else len(read)


@pytest.mark.parametrize("case", CASES, ids=describe)
def test_where_the_gap_is_written_does_not_change_the_result(tmp_path, case):
    kind, gap, run = case
    reference, read, cigars = homopolymer(kind, gap, run)
    data = placements(tmp_path, "ambiguous", reference, read, cigars)

    run_cmuts(data, tmp_path / "banded.h5", band=gap, **WEIGHTED)
    written = written_rows(tmp_path / "banded.h5", len(cigars))

    # Rows of zeros would agree with one another, so the comparison says nothing
    # unless every placement counted the read it was given.
    covered = np.nansum(written["coverage"], axis=1)
    assert np.allclose(covered, scored(kind, gap, read), rtol=0.05, atol=0)

    for name, rows in written.items():
        for cigar, row in zip(cigars, rows):
            assert np.allclose(row, rows[0], atol=TOLERANCE, rtol=0), \
                f"{name} at {cigar} differs from {cigars[0]}"


@pytest.mark.parametrize("case", CONTROLS, ids=describe)
def test_a_band_narrower_than_the_gap_leaves_the_placements_apart(tmp_path, case):
    """What the band is doing, stated the other way around.

    Without it the answer is whatever the aligner chose, and the rows stand as
    far apart as the event is large. The mutations say so and the coverage does
    not: an inserted base covers nothing wherever it is placed, so that array
    can agree across placements even where nothing has been marginalized.
    """
    kind, gap, run = case
    reference, read, cigars = homopolymer(kind, gap, run)
    data = placements(tmp_path, "ambiguous", reference, read, cigars)

    run_cmuts(data, tmp_path / "narrow.h5", band=gap - 1, **WEIGHTED)

    narrow = written_rows(tmp_path / "narrow.h5", len(cigars))

    assert spread(narrow["reactivity"]) > DIVIDED
