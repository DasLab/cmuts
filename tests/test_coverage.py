"""Coverage at a band of zero agrees with samtools depth.

A band of zero pins the marginal to the CIGAR as written, so each kept read
covers exactly the positions its aligned bases pair. samtools depth counts the
same positions, once the reads cmuts hmm always rejects are removed, along
with the reads whose CIGAR the model cannot follow as written.
"""

import numpy as np

from alignments import drop_adjacent_indels
from filters import UNFILTERED
from oracle import rows_by_name, samtools_depth
from outputs import COVERAGE, field_of
from programs import run_cmuts

# How far a stored coverage may stand from the integer depth. Each read's
# contribution is a posterior of one, accumulated in doubles and stored as a
# 32-bit float, so the gap is storage rounding and not model spread.
RELATIVE = 1e-6
ABSOLUTE = 1e-4


def expected_row(depths: dict, columns: int) -> np.ndarray:
    """Builds the coverage row a reference should hold: its depth at every
    covered position, and the fill of zero everywhere else."""
    row = np.zeros(columns)

    for position, depth in depths.items():
        row[position] = depth

    return row


def test_coverage_at_band_zero_matches_samtools_depth(data, falsifiable, tmp_path):
    pinned = drop_adjacent_indels(data, tmp_path)
    depths = samtools_depth(pinned, tmp_path)

    falsifiable(any(depths.values()))

    output = run_cmuts(pinned, tmp_path / "out.h5", band=0, **UNFILTERED)
    coverage = field_of(output, COVERAGE)

    for name, row in rows_by_name(pinned.fasta).items():
        expected = expected_row(depths[name], coverage.shape[1])

        assert np.allclose(coverage[row], expected,
                           rtol=RELATIVE, atol=ABSOLUTE), name
