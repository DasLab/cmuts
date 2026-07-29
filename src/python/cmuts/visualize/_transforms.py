"""Pure-numpy data preparation for the plotly plot builders.

Every plot is "compute the array(s) to draw, then draw them". The compute half
lives here -- numpy in, numpy out, no plotly imports -- so the drawing code stays
thin and the transforms are unit-testable on their own.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
import numpy.typing as npt

# Correlation matrices use a symmetric-log colour scale around zero.
CORRELATION_LINTHRESH = 1e-3


def termination_density(term: npt.ArrayLike) -> np.ndarray:
    """Mean per-position termination density: row-normalize each reference's
    termination counts, then average over references.
    """
    term = np.asarray(term, dtype=float).copy()
    quot = term.sum(axis=1)[:, None]
    term = np.divide(term, quot, where=(quot > 0), out=term)
    return np.mean(term, axis=0)


def reads_per_block(reads: npt.ArrayLike, nblocks: int = 100) -> np.ndarray:
    """Mean reads per reference, references grouped into ``nblocks`` equal-size
    blocks in FASTA order (no sorting; preserves any spatial trend).
    """
    reads = np.asarray(reads)
    nblocks = max(min(nblocks, reads.shape[0]), 1)
    block = max(reads.shape[0] // nblocks, 1)
    reads = reads[: nblocks * block]
    return reads.reshape(nblocks, block).mean(1)


def coverage_fraction(coverage: npt.ArrayLike, reads: npt.ArrayLike) -> np.ndarray:
    """Per-position coverage as a fraction of the mean read depth."""
    return np.asarray(coverage) / np.asarray(reads).mean()


def read_log_depth(reads: npt.ArrayLike) -> np.ndarray:
    """log10 read depth per reference; zero-read references map to -1."""
    reads = np.asarray(reads)
    with np.errstate(divide="ignore"):
        return np.where(reads == 0, -1, np.log10(reads))


def read_histogram(
    reads: npt.ArrayLike, bins: int = 100
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Histogram of log10 read depth.

    Returns ``(bin_centers, counts, normalized)``, where ``normalized`` rescales
    the counts to [0, 1] for per-bar colouring (zeros when all counts are equal).
    Each backend draws and colours from these; only the binning is shared.
    """
    counts, edges = np.histogram(read_log_depth(reads), bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    cmin, cmax = float(counts.min()), float(counts.max())
    if cmax == cmin:
        normalized = np.zeros_like(counts, dtype=float)
    else:
        normalized = (counts - cmin) / (cmax - cmin)
    return centers, counts, normalized


def symlog(x: npt.ArrayLike, linthresh: float = CORRELATION_LINTHRESH) -> np.ndarray:
    """Symmetric log: ``sign(x) * log10(1 + |x| / linthresh)``."""
    x = np.asarray(x)
    return np.sign(x) * np.log10(1.0 + np.abs(x) / linthresh)


# Colour-scale bounds for the log-scaled matrix plots (plotly applies log10 by
# hand). A NamedTuple so callers can either unpack ``vlow, vhigh`` or use ``.vlow``.
class LogBounds(NamedTuple):
    vlow: float
    vhigh: float


def mi_bounds(values: npt.ArrayLike) -> LogBounds:
    """Colour bounds for a mutual-information matrix: the 95th/99th percentile
    of the non-NaN entries, each floored at 1e-6.
    """
    values = np.asarray(values)
    mask = ~np.isnan(values)
    return LogBounds(
        max(float(np.percentile(values[mask], 95)), 1e-6),
        max(float(np.percentile(values[mask], 99)), 1e-6),
    )


def pairwise_coverage_bounds(values: npt.ArrayLike) -> LogBounds:
    """Colour bounds for a pairwise-coverage matrix: [max(min, 1e-4), 1]."""
    return LogBounds(max(float(np.asarray(values).min()), 1e-4), 1.0)
