"""Marginalizing over the band, where the aligner's choice of CIGAR is arbitrary.

A gap inside a homopolymer can be written after any base of the run, and every
CIGAR that writes it describes the same alignment. The arrays must therefore
come out the same for every one of them, which holds once the band is as wide
as the gap, at any length of run.
"""

import h5py
import numpy as np
import pytest

from support import datasets_of, placements, run_cmuts

KINDS = ("D", "I")
GAPS = (1, 2, 3)
RUNS = (10, 30)

# Flanks pin both ends of the read: sliding it by a base costs two mismatches
# against the run, more than the slide saves. They differ from each other so
# that neither end of the read can match the wrong one.
LEFT = "GC"
RIGHT = "TG"

# How far two rows may stand apart and still be the same answer. A band is drawn
# around each row's own CIGAR, so placements admit slightly different paths at
# the edges. The difference is real rather than rounding: around 6e-5, three
# orders above the resolution the counts are stored at, shrinking as the band
# widens.
TOLERANCE = 1e-3

# How far apart a band too narrow for its gap must leave them. Well inside the
# smallest departure measured, which is a third of a mutation.
DIVIDED = 0.1

# Insertions carry no weight by default, which would leave every insertion case
# comparing arrays of zeros. One read to a reference leaves the evidence below a
# whole observation, which would leave every rate NaN.
WEIGHTED = dict(insertion_weight=1, min_depth=0)

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
    """Every dataset a row can be read from, by name, checked for one row per
    placement.

    The row count is checked here because a comparison of rows to one another
    passes over however many rows there are.
    """
    with h5py.File(output, "r") as handle:
        # A run total has no rows to compare, being one number for the file.
        rows = {k: d[:] for k, d in datasets_of(handle).items() if d.ndim >= 1}

    for name, values in rows.items():
        assert len(values) == placements, \
            f"{name}: {len(values)} rows for {placements} placements"

    return rows


def spread(rows) -> float:
    """How far apart the rows of one dataset stand, over every position."""
    return float(np.nanmax(np.nanmax(rows, axis=0) - np.nanmin(rows, axis=0)))


def scored(kind: str, gap: int, read: str) -> int:
    """Bases of the read that meet a reference position, which excludes an
    inserted base."""
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
            assert np.allclose(row, rows[0], atol=TOLERANCE, rtol=0,
                               equal_nan=True), \
                f"{name} at {cigar} differs from {cigars[0]}"


@pytest.mark.parametrize("case", CONTROLS, ids=describe)
def test_a_band_narrower_than_the_gap_leaves_the_placements_apart(tmp_path, case):
    """Asserted over the reactivity and not the coverage: an inserted base
    covers nothing at any placement, so the coverage can agree across
    placements even when nothing has been marginalized."""
    kind, gap, run = case
    reference, read, cigars = homopolymer(kind, gap, run)
    data = placements(tmp_path, "ambiguous", reference, read, cigars)

    run_cmuts(data, tmp_path / "narrow.h5", band=gap - 1, **WEIGHTED)

    narrow = written_rows(tmp_path / "narrow.h5", len(cigars))

    assert spread(narrow["reactivity"]) > DIVIDED
