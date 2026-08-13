"""What the weights on each kind of difference count towards.

A weight scales what its kind of difference contributes to the mutation total.
Nothing here asserts what any rate comes to, only what holds of the rates under
a setting of the weights.

The insertion weight also enters the denominator, an inserted base occupying no
reference position, so a rate is proportional to the weights only where that one
is zero.
"""

import h5py
import numpy as np
import pytest

from datasets import DATASETS
from support import run_cmuts

# Every kind of difference a read can carry, and what each is worth by default.
WEIGHTS = ("substitution_weight", "deletion_weight", "insertion_weight")

# Powers of two, which scale a stored rate without rounding it, so the
# comparison is exact.
SCALES = [0.5, 0.25]

# The set of positions carrying a rate is held fixed by asking for no evidence:
# what a position needs counts the weighted differences, so a lower weight
# otherwise leaves some of them without a rate to compare.
EVERY_POSITION = dict(min_depth=0)


def rates(output):
    """The reactivity and its error, at the positions carrying one."""
    with h5py.File(output, "r") as handle:
        reactivity = handle["reactivity"][:]
        error = handle["error"][:]

    known = ~np.isnan(reactivity)

    return reactivity[known], error[known]


def reactivity(output):
    """The rate at every position, missing values and all."""
    with h5py.File(output, "r") as handle:
        return handle["reactivity"][:]


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_weighing_every_difference_at_zero_leaves_nothing_counted(datasets, checked,
                                                                  tmp_path, name):
    """A difference worth nothing reaches no total, whatever the alignment says
    it is: the rates come to zero rather than to something small."""
    data = datasets(name)
    output = tmp_path / f"{name}.h5"

    run_cmuts(data, output, **dict.fromkeys(WEIGHTS, 0))

    reactivity, error = rates(output)

    checked(reactivity)

    assert (reactivity == 0).all(), f"a rate of {reactivity.max()} with nothing weighed"
    assert (error == 0).all(), f"an error of {error.max()} with nothing weighed"


@pytest.mark.parametrize("scale", SCALES)
@pytest.mark.parametrize("name", sorted(DATASETS))
def test_scaling_the_weights_scales_every_rate_by_the_same_factor(
    datasets, checked, tmp_path, name, scale,
):
    """Insertions are weighed at zero throughout, they being the kind that
    reaches the denominator as well."""
    data = datasets(name)
    weighed = dict(insertion_weight=0, **EVERY_POSITION)

    full = tmp_path / "full.h5"
    scaled = tmp_path / "scaled.h5"

    run_cmuts(data, full, substitution_weight=1, deletion_weight=1, **weighed)
    run_cmuts(data, scaled, substitution_weight=scale, deletion_weight=scale, **weighed)

    before, after = reactivity(full), reactivity(scaled)
    known = ~np.isnan(before)

    assert np.array_equal(np.isnan(before), np.isnan(after)), \
        "the positions carrying a rate are not the same ones"
    assert np.array_equal(after[known], scale * before[known]), \
        "a rate is not the weights' own multiple of itself"

    # A rate of zero scales to zero whatever the factor, so only the positions
    # carrying one say anything about the weights.
    checked(before[known][before[known] > 0])
