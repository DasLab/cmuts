"""Tools for computing what cmuts norm should divide by, and what dividing by
it leaves.

Both schemes pool the aggregate rate of every input and return one number.
The aggregate of a position is one less the product of the pooled channels'
no-event rates, computed in float32. The norm is computed over the pool in
float64, as the program does, so a caller should allow a tolerance.
"""

from __future__ import annotations

import math

import numpy as np

from outputs import (COVERAGE, ERROR_FIELDS, POOLED_RATE_FIELDS, RATE_FIELDS,
                     field_of)

UBR = "ubr"
OUTLIER = "outlier"

# The rate the ubr norm sits at, as a percentile of the pool.
UBR_PERCENTILE = 90

# The band the outlier norm averages, as fractions of the pool counted from the
# highest value down.
OUTLIER_HIGHEST = 0.02
OUTLIER_LOWEST = 0.10

# The coverage a position needs before its rate joins the pool.
MIN_COVERAGE = 500


def _pooled(name, inputs):
    """Every value of one field, from every input, in one array."""
    return np.concatenate([np.asarray(field_of(path, name)).ravel() for path in inputs])


def _aggregate(inputs):
    """The rate of an event of any kind at each position: one less the product
    of the pooled channels' no-event rates, in float32 as the program computes
    it. NaN in any of them carries through."""
    none = np.float32(1.0)

    for name in POOLED_RATE_FIELDS:
        none = none * (np.float32(1.0) - _pooled(name, inputs).astype(np.float32))

    return np.float32(1.0) - none


def pool(inputs, min_coverage=MIN_COVERAGE):
    """Gives the values the norm comes from, which are the aggregates of the
    positions whose coverage clears the floor."""
    rates = _aggregate(inputs)
    keep = np.isfinite(rates) & (_pooled(COVERAGE, inputs) > min_coverage)

    return rates[keep]


def _rank_from_top(n: int, fraction: float) -> int:
    """The rank a fraction of the way down from the highest value, counting from
    zero. Rounds half away from zero, as the program does."""
    return max(1, math.floor((n * fraction) + 0.5) - 1)


def _ubr_norm(values) -> float:
    if values.size == 0:
        return 1.0

    return float(np.percentile(values.astype(np.float64), UBR_PERCENTILE))


def _outlier_norm(values) -> float:
    if values.size < 2:
        return 1.0

    highest = _rank_from_top(values.size, OUTLIER_HIGHEST)
    lowest = _rank_from_top(values.size, OUTLIER_LOWEST)
    ranked = np.sort(values.astype(np.float64))[::-1]

    return float(np.mean(ranked[highest:lowest + 1]))


def norm(scheme, inputs, min_coverage=MIN_COVERAGE) -> float:
    """What cmuts norm should divide every input by. Where the pool supports no
    norm, the value is one and the rates are left as they are."""
    values = pool(inputs, min_coverage)
    found = _ubr_norm(values) if scheme == UBR else _outlier_norm(values)

    return 1.0 if math.isnan(found) or found <= 0 else found


def normalized(path, name, norm) -> np.ndarray:
    """One field of an input divided by the norm, rounded as the program rounds
    it: each value widened to a double for the division and narrowed back."""
    values = np.asarray(field_of(path, name)).astype(np.float64)

    return (values / norm).astype(np.float32)


def expected(path, name, norm) -> np.ndarray:
    """What one field of an output should hold. Every rate and every error
    carries the norm; everything else is copied."""
    if name not in RATE_FIELDS + ERROR_FIELDS:
        return np.asarray(field_of(path, name))

    return normalized(path, name, norm)
