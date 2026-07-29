"""Normalization schemes for reactivity profiles.

A scheme defines only a *formula*: how to turn a pool of reactivity values into a
scale factor. The *granularity* of that pool is chosen separately, by the two
axes on :class:`~cmuts.internal.Opts` (``per_experiment``, ``per_reference``), so
any scheme composes with any granularity. Schemes are small classes registered in
a ``name -> instance`` table; both the CLI ``--norm`` choices and the dispatch
read from it, so adding a scheme is a single localized change (subclass
:class:`Scheme`, decorate with :func:`register`, implement :meth:`block_factor`).

Built-in schemes:
    raw     : no normalization
    ubr     : upper-bound reference -- 90th percentile of high-coverage positions
    outlier : 2-8% outlier-based normalization
    sm-dms  : ShapeMapper2-style per-nucleotide DMS normalization (per-base 75th
              percentile)
    sm-shape: ShapeMapper2-style SHAPE boxplot normalization

Granularity (one normalization factor is computed per "pool" of reference
profiles, then divided out):

    per_experiment=False : pool across all experiments (one shared factor keeps
                           experiments on a comparable scale)
    per_reference=False  : pool across all references (one factor for the pool)

So with neither axis set there is a single factor over everything; with both set
each (experiment, reference) is normalized independently (like rf-norm).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..internal import Opts, ProbingData

__all__ = [
    "Scheme",
    "normalization",
    "register",
    "scheme_names",
]


# UBR (upper-bound reference) normalization parameters.
_UBR_READ_CUTOFF = 500  # per-reference read floor for "high-coverage" positions
_UBR_PERCENTILE = 90  # reactivity percentile used as the normalization factor

# 2-8% (outlier) normalization band.
_OUTLIER_LOW = 0.02  # drop the top 2% as outliers
_OUTLIER_HIGH = 0.10  # average down to the top 10%

# ShapeMapper2-style per-nucleotide DMS normalization parameters.
_DMS_PERCENTILE = 75  # per-base reactivity percentile used as the factor
_DMS_MIN_FACTOR = 0.002  # below this, a base's signal is treated as absent (factor -> NaN)
_DMS_BASE_TOKENS = (0, 1, 2, 3)  # A, C, G, U tokens (see internal._BASE_TOKEN)

# ShapeMapper2-style SHAPE boxplot normalization parameters.
_SM_SHAPE_MIN = 20  # minimum finite values to attempt a factor
_SM_SHAPE_SMALL = 100  # below this many values, exclude the top 5% (not 10%) as outliers


# ---------------------------------------------------------------------------
# Scheme base class and registry
# ---------------------------------------------------------------------------


class Scheme(ABC):
    """A normalization scheme: the factor formula for a pool of reference profiles.

    Subclasses set ``name`` and implement :meth:`block_factor`, which receives one
    pool (the reference profiles being normalized together under the chosen
    granularity) and returns a per-position scale factor for it. The granularity
    -- which profiles share a pool -- is handled by :func:`normalization`, so a
    scheme never needs to know about experiments or references.
    """

    name: str
    needs_sequence: bool = False  # True if block_factor reads the base tokens

    @abstractmethod
    def block_factor(
        self, react: np.ndarray, reads: np.ndarray, bases: np.ndarray | None
    ) -> np.ndarray:
        """Per-position scale factor for one pool of ``n`` reference profiles.

        Args:
            react: ``(n, P)`` reactivities for the pooled profiles.
            reads: ``(n,)`` per-profile read count (per-reference coverage).
            bases: ``(n, P)`` base tokens, or ``None`` for sequence-agnostic schemes.

        Returns:
            ``(n, P)`` factor array. Uniform schemes broadcast one scalar; per-base
            schemes vary it by base. NaN entries leave that position unscaled.
        """


_SCHEMES: dict[str, Scheme] = {}


def register(scheme_cls: type[Scheme]) -> type[Scheme]:
    """Class decorator: register an instance of ``scheme_cls`` under its name."""
    instance = scheme_cls()
    _SCHEMES[instance.name] = instance
    return scheme_cls


def scheme_names() -> list[str]:
    """Names of all registered schemes (drives the CLI ``--norm`` choices)."""
    return list(_SCHEMES)


def _resolve(name: str) -> Scheme:
    try:
        return _SCHEMES[name.lower()]
    except KeyError:
        raise ValueError(
            f"Unknown normalization scheme: {name}. Use one of {scheme_names()}."
        ) from None


# ---------------------------------------------------------------------------
# Built-in schemes
# ---------------------------------------------------------------------------


def _uniform(react: np.ndarray, factor: float) -> np.ndarray:
    """A pool-wide scalar broadcast to every position."""
    return np.full(react.shape, factor, dtype=float)


@register
class RawScheme(Scheme):
    """No normalization."""

    name = "raw"

    def block_factor(self, react, reads, bases):
        return np.ones(react.shape, dtype=float)


@register
class UBRScheme(Scheme):
    """Upper-bound reference: the percentile of reactivity over high-coverage
    positions in the pool. Falls back to no normalization when no position clears
    the coverage floor.
    """

    name = "ubr"

    def block_factor(self, react, reads, bases):
        vals = react[reads > _UBR_READ_CUTOFF]
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            return np.ones(react.shape, dtype=float)
        return _uniform(react, float(np.percentile(vals, _UBR_PERCENTILE)))


@register
class OutlierScheme(Scheme):
    """2-8% normalization: from the top 10% of the pool, drop the top 2% and
    divide by the mean of the remaining 8%.
    """

    name = "outlier"

    def block_factor(self, react, reads, bases):
        vals = react[np.isfinite(react)]
        if vals.size == 0:
            return np.ones(react.shape, dtype=float)
        p_low = max(1, round(vals.size * _OUTLIER_LOW) - 1)
        p_high = max(1, round(vals.size * _OUTLIER_HIGH) - 1)
        ranked = np.sort(vals)[::-1]
        return _uniform(react, float(np.mean(ranked[p_low : p_high + 1])))


@register
class ShapeMapperDMSScheme(Scheme):
    """ShapeMapper2-style per-nucleotide DMS normalization.

    Each base type (A, C, G, U) is scaled independently by the 75th percentile of
    that base's reactivities over the pool (matching ShapeMapper2's per-nucleotide
    scheme). A position is divided by its own base's factor. A base whose factor
    falls below ``_DMS_MIN_FACTOR`` is treated as having no usable signal (factor
    NaN -> left unscaled by :meth:`ProbingData.normalize`).

    Requires the per-position sequence, which ``compute_reactivity`` attaches from
    the FASTA when the counts file carries none.
    """

    name = "sm-dms"
    needs_sequence = True

    def block_factor(self, react, reads, bases):
        out = np.ones(react.shape, dtype=float)
        if bases is None:
            return out
        for tok in _DMS_BASE_TOKENS:
            vals = react[(bases == tok) & np.isfinite(react)]
            f = float(np.percentile(vals, _DMS_PERCENTILE)) if vals.size else float("nan")
            out[bases == tok] = f if f >= _DMS_MIN_FACTOR else float("nan")
        return out


@register
class ShapeMapperSHAPEScheme(Scheme):
    """ShapeMapper2-style SHAPE (2A3) boxplot normalization.

    Reproduces ShapeMapper2's box-plot factor: over the pool, drop outliers above
    ``max(1.5*IQR, the value at the 90th percentile)`` -- whichever cuts fewer
    points -- then divide by the mean of the top 10% of what remains. Below 100
    values the 95th-percentile cut is used in place of the 90th. The factor is a
    single scalar shared by every base, matching ShapeMapper's SHAPE scheme.
    Pools with fewer than 20 finite values are left unnormalized.
    """

    name = "sm-shape"

    def block_factor(self, react, reads, bases):
        vals = np.sort(react[np.isfinite(react)])
        n = vals.size
        ten_pct = n // 10
        if n < _SM_SHAPE_MIN or ten_pct == 0:
            return np.ones(react.shape, dtype=float)
        iqr_limit = 1.5 * abs(np.percentile(vals, 25) - np.percentile(vals, 75))
        tail = ten_pct if n >= _SM_SHAPE_SMALL else n // 20
        limit = max(iqr_limit, vals[n - 1 - tail])
        kept = vals[vals < limit]
        if kept.size < ten_pct:
            return np.ones(react.shape, dtype=float)
        return _uniform(react, float(np.mean(kept[-ten_pct:])))


# ---------------------------------------------------------------------------
# Granularity dispatch
# ---------------------------------------------------------------------------


def normalization(experiments: list[ProbingData], opts: Opts) -> list[np.ndarray]:
    """Per-experiment normalization factors under ``opts.norm`` and the two
    granularity axes (``opts.per_experiment``, ``opts.per_reference``).

    Each experiment's combined reactivity is stacked into ``(E, R, P)``; the
    scheme's factor is then computed once per pool, where a pool spans all
    experiments unless ``per_experiment`` and all references unless
    ``per_reference``. Returns one ``(R, P)`` factor array per experiment, in
    input order, to divide that experiment's reactivity by.
    """
    scheme = _resolve(opts.norm)
    react = np.stack([np.asarray(e.reactivity, dtype=float) for e in experiments])  # (E, R, P)
    reads = np.stack([np.asarray(e.reads, dtype=float) for e in experiments])  # (E, R)
    e_count, r_count, p_count = react.shape

    bases = None
    if scheme.needs_sequence and experiments[0].sequences is not None:
        bases = np.asarray(experiments[0].sequences)  # (R, P), shared across experiments

    exp_pools = [[e] for e in range(e_count)] if opts.per_experiment else [list(range(e_count))]
    ref_pools = [[r] for r in range(r_count)] if opts.per_reference else [list(range(r_count))]

    factors = np.ones((e_count, r_count, p_count), dtype=float)
    for eg in exp_pools:
        for rg in ref_pools:
            idx: tuple[np.ndarray, ...] = np.ix_(eg, rg)
            ne, nr = len(eg), len(rg)
            n = ne * nr
            pool_react = react[idx].reshape(n, p_count)
            pool_reads = reads[idx].reshape(n)
            pool_bases = None
            if bases is not None:
                pool_bases = np.broadcast_to(bases[rg], (ne, nr, p_count)).reshape(n, p_count)
            block = scheme.block_factor(pool_react, pool_reads, pool_bases)
            factors[idx] = block.reshape(ne, nr, p_count)

    return [factors[e] for e in range(e_count)]
