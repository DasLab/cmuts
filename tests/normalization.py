"""Tools for computing the scale cmuts-norm should divide by, and what dividing
by it leaves.

Both schemes pool the rates of every input and return one number. The pool is
gathered in float32 and the scale computed over it in float64, as the program
does, so a caller should allow a tolerance.
"""

from __future__ import annotations

import math

import numpy as np

from outputs import COVERAGE, ERROR, REACTIVITY, field_of

UBR = "ubr"
OUTLIER = "outlier"

# The rate the ubr scale sits at, as a percentile of the pool.
UBR_PERCENTILE = 90

# The band the outlier scale averages, as fractions of the pool counted from the
# highest value down.
OUTLIER_HIGHEST = 0.02
OUTLIER_LOWEST = 0.10

# The coverage a position needs before its rate joins the ubr pool.
MIN_COVERAGE = 500


def _pooled(name, inputs):
    """Every value of one field, from every input, in one array."""
    return np.concatenate([np.asarray(field_of(path, name)).ravel() for path in inputs])


def pool(scheme, inputs, min_coverage=MIN_COVERAGE):
    """The values the scale is taken from. ubr keeps only the positions whose
    coverage clears the floor; outlier keeps every rate there is."""
    rates = _pooled(REACTIVITY, inputs)
    keep = np.isfinite(rates)

    if scheme == UBR:
        keep &= _pooled(COVERAGE, inputs) > min_coverage

    return rates[keep]


def _rank_from_top(n: int, fraction: float) -> int:
    """The rank a fraction of the way down from the highest value, counting from
    zero. Rounds half away from zero, as the program does."""
    return max(1, math.floor((n * fraction) + 0.5) - 1)


def _ubr_factor(values) -> float:
    if values.size == 0:
        return 1.0

    return float(np.percentile(values.astype(np.float64), UBR_PERCENTILE))


def _outlier_factor(values) -> float:
    if values.size < 2:
        return 1.0

    highest = _rank_from_top(values.size, OUTLIER_HIGHEST)
    lowest = _rank_from_top(values.size, OUTLIER_LOWEST)
    ranked = np.sort(values.astype(np.float64))[::-1]

    return float(np.mean(ranked[highest:lowest + 1]))


def factor(scheme, inputs, min_coverage=MIN_COVERAGE) -> float:
    """The scale cmuts-norm should divide every input by. A scale that says
    nothing leaves the rates as they are."""
    values = pool(scheme, inputs, min_coverage)
    found = _ubr_factor(values) if scheme == UBR else _outlier_factor(values)

    return 1.0 if math.isnan(found) or found <= 0 else found


def scaled(path, name, scale) -> np.ndarray:
    """One field of an input divided by the scale, rounded as the program rounds
    it: each value widened to a double for the division and narrowed back."""
    values = np.asarray(field_of(path, name)).astype(np.float64)

    return (values / scale).astype(np.float32)


def clipped(values, below=None, above=None) -> np.ndarray:
    """The values held within the bounds, leaving NaN as it is."""
    held = values

    if below is not None:
        held = np.maximum(held, np.float32(below))
    if above is not None:
        held = np.minimum(held, np.float32(above))

    return held


def expected(path, name, scale, below=None, above=None) -> np.ndarray:
    """What one field of an output should hold. Only the rates and their error
    carry the scale, and only the rates are clipped."""
    if name not in (REACTIVITY, ERROR):
        return np.asarray(field_of(path, name))

    values = scaled(path, name, scale)

    return clipped(values, below, above) if name == REACTIVITY else values
