"""Randomly shaped data, checked the same way as the fixed shapes.

This is where a case neither the fixtures nor the matrix thought to cover turns
up. On a failure hypothesis reports the parameters that produced it and shrinks
them towards the smallest that still fails, so the case arrives ready to
investigate rather than as a large file to bisect by hand.
"""

import os

from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from support import generate, run_cmuts, samtools_kept

EXAMPLES = int(os.environ.get("FUZZ_EXAMPLES", "25"))


def ranges(low, high):
    """A spec of the form LOW:HIGH, with its endpoints the right way round."""
    return st.tuples(st.integers(low, high), st.integers(low, high)).map(
        lambda pair: f"{min(pair)}:{max(pair)}"
    )


def fractions(low, high):
    return st.floats(min_value=low, max_value=high, allow_nan=False, allow_infinity=False)


@st.composite
def shapes(draw):
    """Everything the generator can be asked for."""
    return dict(
        seed=draw(st.integers(0, 2**31)),
        references=draw(st.integers(1, 50)),
        ref_length=draw(ranges(50, 2000)),
        covered=draw(fractions(0.05, 1.0)),
        reads_per_ref=draw(ranges(0, 25)),
        read_length=draw(ranges(20, 500)),
        mapq=draw(ranges(0, 60)),
        reverse=draw(fractions(0.0, 1.0)),
        unmapped=draw(ranges(0, 20)),
        mismatch_rate=draw(fractions(0.0, 0.1)),
        insertions=draw(ranges(0, 2)),
        insertion_length=draw(ranges(1, 150)),
        deletions=draw(ranges(0, 2)),
        deletion_length=draw(ranges(1, 20)),
        soft_clips=draw(ranges(0, 2)),
        soft_clip_length=draw(ranges(1, 60)),
    )


@st.composite
def criteria(draw):
    """A filter, with the upper bound either absent or above the lower one."""
    lower = draw(st.integers(0, 400))

    return dict(
        min_mapq=draw(st.integers(0, 61)),
        strand=draw(st.sampled_from(["both", "forward", "reverse"])),
        min_length=lower,
        max_length=draw(st.one_of(st.just(0), st.integers(lower, lower + 600))),
    )


@settings(
    max_examples=EXAMPLES,
    deadline=None,                                  # every example runs subprocesses
    suppress_health_check=[HealthCheck.too_slow],
)
@given(shape=shapes(), filters=criteria())
def test_random_data_matches_samtools(tmp_path_factory, shape, filters):
    work = tmp_path_factory.mktemp("fuzz")
    data = generate(work, "fuzz", **shape)
    summary = run_cmuts(data, work / "out.h5", **filters)

    assert summary.kept == samtools_kept(data, **filters), "surviving reads"
    assert summary.kept + summary.rejected == data.mapped, "reads accounted for"
    assert summary.unmapped == data.unmapped, "unmapped reads"
    assert summary.rows == data.touched, "references written"
