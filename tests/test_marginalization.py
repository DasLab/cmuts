"""Marginalizing over the band, where the aligner's choice of CIGAR is
arbitrary.

A gap inside a homopolymer can be written after any base of the run, and every
CIGAR that writes it describes the same alignment. The arrays must therefore
come out the same for every one of them, which holds once the band is as wide
as the gap, at any length of run.
"""

from dataclasses import dataclass

import numpy as np
import pytest

from alignments import build_placements
from outputs import COVERAGE, REACTIVITY, ROW_FIELDS, fields_of
from programs import run_cmuts

DELETION = "D"
INSERTION = "I"

GAPS = (1, 2, 3)
RUNS = (10, 30)

# Flanks that pin both ends of the read: sliding it by a base costs two
# mismatches against the run, more than the slide saves. The two differ so that
# neither end of the read can match the wrong flank.
LEFT = "GC"
RIGHT = "TG"

# How far two rows may stand apart and still be the same answer. A band is
# drawn around each row's own CIGAR, so placements admit slightly different
# paths at the edges. The difference is real and not rounding: around 6e-5,
# three orders above the resolution the counts are stored at, and it shrinks as
# the band widens.
TOLERANCE = 1e-3

# How far apart a band too narrow for its gap must leave two rows. Well inside
# the smallest departure measured, which is a third of a mutation.
DIVIDED = 0.1

# Insertions have a weight of zero by default, which would leave every
# insertion case comparing arrays of zeros. A single read leaves the evidence
# below one whole observation, which would leave every rate NaN.
WEIGHTED = dict(insertion_weight=1, min_depth=0)

# How far a placement's coverage may stand from the count of bases the read
# scores. The band allows paths reaching a little further than the CIGAR alone.
COVERED_TOLERANCE = 0.05


@dataclass(frozen=True)
class Ambiguity:
    """A gap of one kind and size, inside a homopolymer run of one length."""

    kind: str
    gap: int
    run: int

    def __str__(self) -> str:
        return f"{self.gap}{self.kind}-in-{self.run}"


CASES = [Ambiguity(kind, gap, run)
         for kind in (DELETION, INSERTION) for gap in GAPS for run in RUNS]

# One run length is enough to show that a narrow band leaves the placements
# apart.
CONTROLS = [case for case in CASES if case.run == RUNS[0]]


def homopolymer(case: Ambiguity):
    """Builds the reference, the read, and every CIGAR that places the gap
    inside the run.

    A deletion takes the gap out of the read and may follow any of the bases it
    leaves behind. An insertion adds the gap, and so has one placement more
    than the run has bases.
    """
    reference = LEFT + "A" * case.run + RIGHT
    added = case.run + case.gap if case.kind == INSERTION else case.run - case.gap
    read = LEFT + "A" * added + RIGHT

    spots = range(case.run + 1) if case.kind == INSERTION else range(case.run - case.gap + 1)
    cigars = []

    for taken in spots:
        left = len(LEFT) + taken
        right = len(read) - left - (case.gap if case.kind == INSERTION else 0)
        cigars.append(f"{left}M{case.gap}{case.kind}{right}M")

    return reference, read, cigars


def read_rows(output, placements: int) -> dict:
    """Reads every field with a row per reference, having checked that there is
    one row per placement.

    The row count is checked here because a comparison of rows to one another
    passes over however many rows there are.
    """
    rows = fields_of(output, ROW_FIELDS)

    for name, values in rows.items():
        assert len(values) == placements, name

    return rows


def row_spread(rows) -> float:
    """Returns how far apart the rows of one field stand, over every
    position."""
    return float(np.nanmax(np.nanmax(rows, axis=0) - np.nanmin(rows, axis=0)))


def scored_bases(case: Ambiguity, read: str) -> int:
    """Returns how many bases of the read meet a reference position, which
    excludes an inserted base."""
    return len(read) - case.gap if case.kind == INSERTION else len(read)


@pytest.mark.parametrize("case", CASES, ids=str)
def test_where_the_gap_is_written_does_not_change_the_result(tmp_path, case):
    reference, read, cigars = homopolymer(case)
    data = build_placements(tmp_path, "ambiguous", reference, read, cigars)

    output = tmp_path / "banded.h5"
    run_cmuts(data, output, band=case.gap, **WEIGHTED)
    written = read_rows(output, len(cigars))

    # Rows of zeros would agree with one another, so the comparison means
    # something only once every placement has counted the read it was given.
    covered = np.nansum(written[COVERAGE], axis=1)
    assert np.allclose(covered, scored_bases(case, read),
                       rtol=COVERED_TOLERANCE, atol=0)

    for name, rows in written.items():
        for cigar, row in zip(cigars, rows):
            assert np.allclose(row, rows[0], atol=TOLERANCE, rtol=0,
                               equal_nan=True), f"{name}: {cigar}"


@pytest.mark.parametrize("case", CONTROLS, ids=str)
def test_a_band_narrower_than_the_gap_leaves_the_placements_apart(tmp_path, case):
    """Asserts over the reactivity: an inserted base covers no reference
    position at any placement, so the coverage can agree across placements even
    where no read has been marginalized."""
    reference, read, cigars = homopolymer(case)
    data = build_placements(tmp_path, "ambiguous", reference, read, cigars)

    output = tmp_path / "narrow.h5"
    run_cmuts(data, output, band=case.gap - 1, **WEIGHTED)

    narrow = read_rows(output, len(cigars))

    assert row_spread(narrow[REACTIVITY]) > DIVIDED
