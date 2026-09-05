"""Tools for computing the scale cmuts norm should divide by, and what dividing
by it leaves.

Both schemes pool the aggregate rate of every input and return one number.
The aggregate of a position is one less the product of the three no-event
rates, computed in float32, and the scale is computed over the pool in float64, as the program
does, so a caller should allow a tolerance.
"""

from __future__ import annotations

import math

import numpy as np

from outputs import COVERAGE, ERROR_FIELDS, RATE_FIELDS, field_of

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


def _aggregate(inputs):
    """The rate of an event of any kind at each position: one less the product
    of the three no-event rates, in float32 as the program computes it. NaN in
    any channel carries through."""
    none = np.float32(1.0)

    for name in RATE_FIELDS:
        none = none * (np.float32(1.0) - _pooled(name, inputs).astype(np.float32))

    return np.float32(1.0) - none


def pool(scheme, inputs, min_coverage=MIN_COVERAGE):
    """The values the scale is taken from. ubr keeps only the positions whose
    coverage clears the floor; outlier keeps every aggregate there is."""
    rates = _aggregate(inputs)
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
    """The scale cmuts norm should divide every input by. Where the pool supports
    no scale, the factor is one and the rates are left as they are."""
    values = pool(scheme, inputs, min_coverage)
    found = _ubr_factor(values) if scheme == UBR else _outlier_factor(values)

    return 1.0 if math.isnan(found) or found <= 0 else found


def scaled(path, name, scale) -> np.ndarray:
    """One field of an input divided by the scale, rounded as the program rounds
    it: each value widened to a double for the division and narrowed back."""
    values = np.asarray(field_of(path, name)).astype(np.float64)

    return (values / scale).astype(np.float32)


def expected(path, name, scale) -> np.ndarray:
    """What one field of an output should hold. Every rate and every error
    carries the scale; everything else is copied."""
    if name not in RATE_FIELDS + ERROR_FIELDS:
        return np.asarray(field_of(path, name))

    return scaled(path, name, scale)
