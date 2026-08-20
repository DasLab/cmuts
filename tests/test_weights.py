"""What the weights on each kind of difference count towards.

A weight scales what its kind of difference contributes to the mutation total.
Each test asserts what holds of the rates under a setting of the weights.

The insertion weight also enters the denominator, an inserted base occupying no
reference position, so a rate is proportional to the weights only where the
insertion weight is zero.
"""

import numpy as np
import pytest

from outputs import ERROR, REACTIVITY, field_of
from programs import run_cmuts

# Every kind of difference a read can carry.
WEIGHTS = ("substitution_weight", "deletion_weight", "insertion_weight")

# Powers of two, which scale a stored rate without rounding it, so the
# comparison is exact.
SCALES = [0.5, 0.25]

# Keeps the set of positions with a rate fixed. The depth requirement counts
# weighted differences, so a lower weight would otherwise leave some positions
# without a rate to compare.
EVERY_POSITION = dict(min_depth=0)


def known_rates(output):
    """Returns the reactivity and its error at the positions where one is defined."""
    reactivity = field_of(output, REACTIVITY)
    error = field_of(output, ERROR)
    known = ~np.isnan(reactivity)

    return reactivity[known], error[known]


def test_weighing_every_difference_at_zero_leaves_the_rates_at_zero(data, falsifiable,
                                                                  tmp_path):
    """A difference of any kind with a weight of zero reaches no total: the
    rates come to zero and not to something small."""
    output = tmp_path / "out.h5"

    run_cmuts(data, output, **dict.fromkeys(WEIGHTS, 0))

    reactivity, error = known_rates(output)

    falsifiable(reactivity.size > 0)

    assert (reactivity == 0).all()
    assert (error == 0).all()


@pytest.mark.parametrize("scale", SCALES)
def test_scaling_the_weights_scales_every_rate_by_the_same_factor(data, falsifiable,
                                                                  tmp_path, scale):
    """Insertions are weighed at zero throughout, being the kind that reaches
    the denominator as well."""
    weighed = dict(insertion_weight=0, **EVERY_POSITION)

    full = tmp_path / "full.h5"
    scaled = tmp_path / "scaled.h5"

    run_cmuts(data, full, substitution_weight=1, deletion_weight=1, **weighed)
    run_cmuts(data, scaled, substitution_weight=scale, deletion_weight=scale, **weighed)

    before, after = field_of(full, REACTIVITY), field_of(scaled, REACTIVITY)
    known = ~np.isnan(before)

    # A rate of zero scales to zero at any factor, so only a rate above zero
    # can detect the scaling.
    falsifiable((before[known] > 0).any())

    assert np.array_equal(np.isnan(before), np.isnan(after))
    assert np.array_equal(after[known], scale * before[known])
