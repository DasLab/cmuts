"""The terminations of a reference sum to the reads counted on it.

Every counted read stops on exactly one base. The rate at a base is that
base's share of the coverage. The reads which stop at a base are therefore
the rate there times the coverage there.

A band of zero pins each read's stop to the 5'-most base of its placed span.
That base always lies inside the reference. A wider band may place the stop
beyond the 5' end, where nothing is accumulated.

The depth floor is lowered so that every position reports a rate. A position
of exactly one read can otherwise fall a rounding step below the floor, which
would drop that read from the sum.
"""

import numpy as np

from alignments import drop_adjacent_indels
from filters import UNFILTERED
from outputs import PRIMARY_COUNTED, COVERAGE, TERMINATION_RATE, field_of
from programs import run_cmuts

# How far a summed posterior may stand from the whole number of reads. Each
# read contributes a posterior of one, accumulated in doubles and stored as a
# 32-bit float, so the gap is storage rounding.
RELATIVE = 1e-6
ABSOLUTE = 1e-4


def test_terminations_at_band_zero_sum_to_the_reads_counted(data, falsifiable, tmp_path):
    pinned = drop_adjacent_indels(data, tmp_path)

    output = run_cmuts(pinned, tmp_path / "out.h5", band=0, min_depth=0, **UNFILTERED)
    rate = field_of(output, TERMINATION_RATE)
    coverage = field_of(output, COVERAGE)
    counted = field_of(output, PRIMARY_COUNTED)

    falsifiable(np.any(np.asarray(counted) > 0))

    stopped = np.nansum(np.asarray(rate) * np.asarray(coverage), axis=1)

    assert np.allclose(stopped, counted, rtol=RELATIVE, atol=ABSOLUTE)
