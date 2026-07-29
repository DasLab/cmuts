"""
Internal implementation of reactivity computation and normalization.

This module contains the core data structures and algorithms for transforming
mutation count data into normalized reactivity profiles. It handles both
traditional 1D MaP-seq analysis and pairwise 2D correlation analysis.

Data Structures:
    ProbingData: Primary container holding reactivity, error, SNR, and pairwise metrics.
    Opts: Configuration dataclass with normalization options (cutoff, blank regions, etc.).
    DataGroups: Paths to HDF5 groups containing count data.
    Datasets: String constants for HDF5 dataset names.

Key Functions:
    compute_reactivity: Main entry point for reactivity calculation.
    _reactivity: Compute raw reactivity as modifications/coverage.
    _error: Compute standard error of Bernoulli proportion.
    _snr: Compute signal-to-noise ratio.
    _correlation: Compute pairwise mutation correlations with significance testing.
    _mutual_information: Compute normalized mutual information matrix.

Data Flow:
    HDF5 counts -> _accumulate_arrays -> _data_from_counts -> ProbingData
    -> _merge_conditions (if nomod) -> normalization -> normalize -> clip -> output
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Union

import dask.array as da
import h5py
import numpy as np
from scipy.stats import t as student

from .normalize.schemes import normalization

# Core datatypes and dataclasses


ProbingDatum = Union[np.ndarray, da.Array]


class DataGroups:
    def __init__(
        self: DataGroups,
        groups: Union[list[str], None],
    ) -> None:
        if groups is None:
            groups = []

        oneD = "/counts-1d"
        twoD = "/counts-2d"

        self.oneD = [group + oneD for group in groups]
        self.twoD = [group + twoD for group in groups]


@dataclass
class Opts:
    mod: DataGroups
    nomod: DataGroups
    cutoff: int
    ins: bool
    dels: bool
    norm: str  # registered scheme name (see normalize.scheme_names)
    blank: tuple[int, int]
    clip: tuple[float | None, float | None]  # (clip_below, clip_above); None = no clip
    sig: float
    # Normalization granularity. A factor is computed per pool of reference
    # profiles; per_experiment pools each experiment separately (else all
    # experiments share a factor), per_reference pools each reference separately
    # (else all references share a factor).
    per_experiment: bool = False
    per_reference: bool = False


@dataclass
class Experiment:
    """A named experiment with one or more replicate datasets.

    Attributes:
        name: Experiment name, used as the output HDF5 group when saving.
        mod: HDF5 paths for the modified condition (replicates are summed).
        nomod: HDF5 paths for the unmodified control (optional).
    """

    name: str
    mod: list[str]
    nomod: Union[list[str], None] = None


# Input array indices


IX_DEL = 4
IX_INS = 5
IX_TERM = 6


# Datasets in the output HDF5 file


@dataclass
class Datasets:
    REACTIVITY = "reactivity"
    READS = "reads"
    ERROR = "error"
    SNR = "SNR"
    PAIRWISE_SNR = "pairwise-snr"
    NORM = "norm"
    ROI = "roi-mask"
    HEATMAP = "heatmap"
    COVERAGE = "coverage"
    TERMINATIONS = "terminations"
    MI = "mutual-information"
    COV = "covariance"
    SEQUENCE = "sequence"
    NAME = "name"
    # Shared, experiment-independent datasets live under this top-level group
    # (SEQUENCE, NAME, ...); keeping them grouped leaves room for future metadata
    # and keeps them from colliding with experiment groups.
    META = "meta"


# Misc utils


def _has_arrays(
    file: h5py.File,
    groups: list[str],
) -> bool:
    if not groups:
        return False

    for group in groups:
        if group not in file:
            return False

    return True


def _accumulate_arrays(
    file: h5py.File,
    groups: list[str],
) -> da.Array:
    arr: da.Array = da.from_array(_get_dataset(file, groups[0]), chunks="auto")

    for group in groups[1:]:
        arr += da.from_array(_get_dataset(file, group), chunks="auto")

    return arr


def _get_dataset(file: h5py.File, path: str) -> h5py.Dataset:
    try:
        return file[path]
    except KeyError:
        available = list(file.keys())
        raise KeyError(
            f"Dataset '{path}' not found in '{file.filename}'. Available groups: {available}"
        ) from None


#
# Traditional MaP-seq
#


def _reactivity(
    modifications: da.Array,
    coverage: da.Array,
    mask: da.Array,
) -> da.Array:
    reactivity = da.divide(
        modifications,
        coverage,
        where=mask,
        out=da.ones_like(modifications) * np.nan,
    )
    return da.clip(reactivity, 0, 1)


def _error(
    reactivity: da.Array,
    coverage: da.Array,
    mask: da.Array,
) -> da.Array:
    return da.divide(
        da.sqrt(reactivity * (1 - reactivity)),
        da.sqrt(coverage),
        where=mask,
        out=da.ones_like(reactivity) * np.nan,
    )


def _snr(
    reactivity: ProbingDatum,
    error: ProbingDatum,
) -> da.Array:
    return da.divide(
        reactivity,
        error,
        where=(error > 0),
        out=np.zeros_like(reactivity),
    ).mean(-1)


def _roi(
    seqs: int,
    seqlen: int,
    opts: Opts,
    lengths: da.Array,
) -> da.Array:
    low = opts.blank[0]
    high = opts.blank[1]

    roi = da.ones((seqs, seqlen), dtype=bool)
    if low:
        roi[:, :low] = False
    if high:
        roi[:, -high:] = False

    length_mask = da.arange(seqlen) < lengths[:, None]
    return roi & length_mask


#
# Pairwise MaP-seq
#


def _validate_pairwise_shape(data: da.Array) -> bool:
    return (
        data.ndim == 5
        and data.shape[1] == data.shape[2]
        and data.shape[3] == 2
        and data.shape[4] == 2
    )


def _probability(pairs: da.Array) -> da.Array:
    if not _validate_pairwise_shape(pairs):
        raise ValueError("The 2D counts have an incorrect shape.")

    prob = da.zeros_like(pairs)
    prob[..., 0, 0] = 1

    counts = pairs.sum((-1, -2))[..., None, None]

    return da.divide(
        pairs,
        counts,
        where=(counts > 0),
        out=prob,
    )


def _conditional(prob: da.Array) -> da.Array:
    joint = prob[..., 1, 1]
    marginal = prob[..., 1, 1] + prob[..., 1, 0]

    return da.divide(joint, marginal, where=(marginal > 0), out=np.zeros_like(joint))


def _diagonal_normalize(
    values: da.Array,
    eps: float = 1e-8,
) -> da.Array:
    diag = da.diagonal(values, axis1=-2, axis2=-1)
    quot = da.sqrt(diag[..., :, None] * diag[..., None, :] + eps)
    out = da.divide(
        values,
        quot,
        where=(quot != 0),
        out=da.zeros_like(values),
    )
    return da.clip(out, -1, 1)


def _covariance(prob: da.Array) -> da.Array:
    ex = prob[..., 1, 1] + prob[..., 1, 0]
    ey = prob[..., 1, 1] + prob[..., 0, 1]
    exy = prob[..., 1, 1]

    return exy - ex * ey


def _sig_test(
    corr: da.Array,
    reads: da.Array,
    sig: float,
    eps: float = 1e-8,
) -> da.Array:
    # Bias correction

    reads = da.maximum(reads - 2, 0)

    # Calculate t-statistic and p-values

    inter = reads / (1 - corr**2 + eps)
    t = corr * da.sqrt(inter)
    p = 2 * (1 - da.from_array(student.cdf(da.abs(t), reads)))

    # Bonferroni correction

    n = corr.shape[-1]
    tests = n * (n - 1) // 2
    p = da.minimum(p * tests, 1.0)

    return p <= sig


def _correlation(
    prob: da.Array,
    reads: da.Array,
    p: Union[float, None],
) -> da.Array:
    cov = _covariance(prob)
    corr = _diagonal_normalize(cov)

    if p is not None:
        ix = _sig_test(corr, reads, p)
        corr = da.where(ix, corr, np.nan)

    return corr


def _mutual_information(prob: da.Array, norm: bool = True) -> da.Array:
    p_00 = prob[..., 0, 0]
    p_01 = prob[..., 0, 1]
    p_10 = prob[..., 1, 0]
    p_11 = prob[..., 1, 1]

    p_x0 = p_00 + p_01
    p_x1 = p_10 + p_11
    p_y0 = p_00 + p_10
    p_y1 = p_01 + p_11

    # p(0,0) * log(p(0,0) / (p(0) * p(0)))
    valid_00 = (p_00 > 0) & (p_x0 > 0) & (p_y0 > 0)
    ratio_00 = da.divide(p_00, p_x0 * p_y0, where=valid_00, out=np.ones_like(p_00))
    term_00 = da.where(valid_00, p_00 * da.log(ratio_00), 0)

    # p(0,1) * log(p(0,1) / (p(0) * p(1)))
    valid_01 = (p_01 > 0) & (p_x0 > 0) & (p_y1 > 0)
    ratio_01 = da.divide(p_01, p_x0 * p_y1, where=valid_01, out=np.ones_like(p_01))
    term_01 = da.where(valid_01, p_01 * da.log(ratio_01), 0)

    # p(1,0) * log(p(1,0) / (p(1) * p(0)))
    valid_10 = (p_10 > 0) & (p_x1 > 0) & (p_y0 > 0)
    ratio_10 = da.divide(p_10, p_x1 * p_y0, where=valid_10, out=np.ones_like(p_10))
    term_10 = da.where(valid_10, p_10 * da.log(ratio_10), 0)

    # p(1,1) * log(p(1,1) / (p(1) * p(1)))
    valid_11 = (p_11 > 0) & (p_x1 > 0) & (p_y1 > 0)
    ratio_11 = da.divide(p_11, p_x1 * p_y1, where=valid_11, out=np.ones_like(p_11))
    term_11 = da.where(valid_11, p_11 * da.log(ratio_11), 0)

    mi = term_00 + term_01 + term_10 + term_11
    if norm:
        return _diagonal_normalize(mi)
    else:
        return mi


def _pairwise_snr(
    prob: da.Array,
    pairs: da.Array,
) -> da.Array:
    """Compute SNR for pairwise joint probability P(i=1, j=1).

    SNR = P(1,1) / SE(P(1,1))
    where SE = sqrt(p * (1-p) / n) for Bernoulli proportion.

    Returns mean SNR across all position pairs for each reference.
    """
    p11 = prob[..., 1, 1]
    n = pairs.sum((-1, -2))

    # Standard error of Bernoulli proportion
    se = da.sqrt(
        da.divide(
            p11 * (1 - p11),
            n,
            where=(n > 0),
            out=da.ones_like(p11),
        )
    )

    # SNR = signal / noise
    snr = da.divide(
        p11,
        se,
        where=(se > 0),
        out=da.zeros_like(p11),
    )

    # Return mean SNR across position pairs (excluding diagonal)
    # Mask out diagonal (i == j)
    L = snr.shape[-1]
    diag_mask = da.eye(L, dtype=bool)
    snr_masked = da.where(diag_mask, 0, snr)
    n_pairs = L * (L - 1)

    return snr_masked.sum(axis=(-1, -2)) / n_pairs


#
# Main normalization function
#


def _pairwise_mask(mask: da.Array) -> da.Array:
    return mask[..., :, None] & mask[..., None, :]


def _data_from_counts(
    counts: da.Array, pairs: Union[da.Array, None], opts: Opts, lengths: da.Array
) -> tuple[Union[da.Array, None], ...]:
    seqs = int(counts.shape[0])
    seqlen = int(counts.shape[1])

    roi = _roi(seqs, seqlen, opts, lengths)

    # Compute the modification heatmap
    # The termination events are normalized separately

    heatmap = da.mean(counts, axis=(0, 1))
    heatmap = heatmap / heatmap.sum()

    # Per-reference coverage: matches + mismatches at each position. Stored per
    # reference (n_ref, L) so the report can plot both the aggregate (mean over
    # references) and a single sequence's coverage.

    coverage = da.sum(counts[..., :IX_DEL], axis=(2, 3))

    # Compute the coverage mask

    mask = coverage >= opts.cutoff
    mask[~roi] = False

    # Compute termination events

    term = da.sum(counts[..., IX_TERM], axis=2)

    # Compute the trace to get the non-modified rates only

    trace = da.trace(counts, axis1=2, axis2=3)

    # Subtract the trace to get the pure modifications

    modifications = coverage - trace

    # Account for insertions and deletions

    if opts.ins:
        modifications += da.sum(counts[..., IX_INS], axis=-1)
    if opts.dels:
        modifications += da.sum(counts[..., IX_DEL], axis=-1)

    # Divide by the coverage to get the un-normalised reactivity

    reactivity = _reactivity(modifications, coverage, mask)

    # SEM of a Bernoulli random variable

    err = _error(reactivity, coverage, mask)

    # Get the read counts and compute all values

    reads = da.max(coverage, axis=(-1,))

    # Signal-to-noise

    snr = _snr(reactivity, err)

    # Pairwise probabilities

    if pairs is not None:
        prob = _probability(pairs)

        mask_2d = _pairwise_mask(mask)

        # Compute pairwise SNR before summing pairs
        pairwise_snr = _pairwise_snr(prob, pairs)

        pairs = pairs.sum((-2, -1))
        assert pairs is not None  # narrowing for mypy
        covariance = _correlation(prob, pairs, opts.sig)
        mi = _mutual_information(prob)

        covariance[~mask_2d] = np.nan
        mi[~mask_2d] = np.nan

    else:
        prob = None
        covariance = None
        mi = None
        pairwise_snr = None

    return (
        reactivity,
        reads,
        err,
        snr,
        mask,
        heatmap,
        coverage,
        term,
        pairs,
        prob,
        covariance,
        mi,
        pairwise_snr,
    )


class ProbingData:
    def __init__(
        self: ProbingData,
        sequences: Union[ProbingDatum, None],
        reactivity: ProbingDatum,
        reads: ProbingDatum,
        error: ProbingDatum,
        snr: ProbingDatum,
        mask: ProbingDatum,
        heatmap: ProbingDatum,
        coverage: ProbingDatum,
        terminations: ProbingDatum,
        pairs: Union[ProbingDatum, None],
        probability: Union[ProbingDatum, None],
        covariance: Union[ProbingDatum, None],
        mi: Union[ProbingDatum, None],
        pairwise_snr: Union[ProbingDatum, None] = None,
    ) -> None:
        self.sequences = sequences
        self.norm: Union[ProbingDatum, None] = None

        # 1D values

        self.reactivity = reactivity
        self.reads = reads
        self.error = error
        self.snr = snr
        self.mask = mask
        self.heatmap = heatmap
        self.coverage = coverage
        self.terminations = terminations

        # 2D values

        self.pairs = pairs
        self.probability = probability
        self.covariance = covariance
        self.mi = mi
        self.pairwise_snr = pairwise_snr

    def size(
        self: ProbingData,
    ) -> int:
        return int(self.reads.shape[0])

    def single(
        self: ProbingData,
    ) -> bool:
        return self.size() == 1

    def __repr__(
        self: ProbingData,
    ) -> str:
        from .output import _print_stats_single

        return _print_stats_single(self)

    def compute(
        self: ProbingData,
    ) -> ProbingData:
        computed = da.compute(
            self.sequences,
            self.reactivity,
            self.reads,
            self.error,
            self.snr,
            self.mask,
            self.heatmap,
            self.coverage,
            self.terminations,
            self.pairs,
            self.probability,
            self.covariance,
            self.mi,
            self.pairwise_snr,
        )
        return ProbingData(*computed)

    def clip(
        self: ProbingData,
        below: float | None,
        above: float | None,
    ) -> None:
        if below is not None:
            self.reactivity = da.clip(self.reactivity, below, None)
        if above is not None:
            self.reactivity = da.clip(self.reactivity, None, above)
        if below is not None or above is not None:
            self.snr = _snr(self.reactivity, self.error)

    def normalize(
        self: ProbingData,
        norm: ProbingDatum,
    ) -> None:
        norm_arr = np.nan_to_num(np.asarray(norm), nan=1)

        if norm_arr.ndim > 0:
            norm_arr = np.where(norm_arr <= 0, 1, norm_arr)
        elif norm_arr <= 0:
            norm_arr = np.ones(1, dtype=norm_arr.dtype)

        self.reactivity /= norm_arr
        self.error /= norm_arr
        self.norm = norm_arr

    def datasets(self: ProbingData) -> dict[str, np.ndarray]:
        """Per-group HDF5 datasets (name -> array).

        The single source of which datasets persist; both :meth:`save` and
        :func:`save_groups` use it, so they cannot drift. Optional fields are
        included only when set. Sequences are file-level and handled separately.
        """
        data: dict[str, np.ndarray] = {
            Datasets.REACTIVITY: np.asarray(self.reactivity),
            Datasets.READS: np.asarray(self.reads),
            Datasets.ERROR: np.asarray(self.error),
            Datasets.SNR: np.asarray(self.snr),
            Datasets.HEATMAP: np.asarray(self.heatmap),
            Datasets.COVERAGE: np.asarray(self.coverage),
            Datasets.TERMINATIONS: np.asarray(self.terminations),
        }
        for key, value in (
            (Datasets.NORM, self.norm),
            (Datasets.MI, self.mi),
            (Datasets.COV, self.covariance),
            (Datasets.PAIRWISE_SNR, self.pairwise_snr),
        ):
            if value is not None:
                data[key] = np.asarray(value)
        return data

    def save(
        self: ProbingData,
        group: str,
        file: str,
    ) -> None:
        data: dict[str, da.Array] = {
            (f"{group}/{key}" if group else key): da.from_array(arr)
            for key, arr in self.datasets().items()
        }
        if self.sequences is not None:
            data[f"{Datasets.META}/{Datasets.SEQUENCE}"] = da.from_array(self.sequences)
        da.to_hdf5(file, data)

    @classmethod
    def load(
        cls,
        group: str,
        file: h5py.File,
    ) -> ProbingData:
        """Load ProbingData from an HDF5 file."""

        def get(name: str) -> Union[np.ndarray, None]:
            key = f"{group}/{name}" if group else name
            return np.array(file[key]) if key in file else None

        # Sequences live under the shared "meta" group; fall back to the old
        # top-level dataset for files written before that move.
        _seq_key = f"{Datasets.META}/{Datasets.SEQUENCE}"
        if _seq_key in file:
            sequences = np.array(file[_seq_key])
        elif Datasets.SEQUENCE in file:
            sequences = np.array(file[Datasets.SEQUENCE])
        else:
            sequences = None
        reactivity = get(Datasets.REACTIVITY)
        reads = get(Datasets.READS)
        error = get(Datasets.ERROR)
        snr = get(Datasets.SNR)
        heatmap = get(Datasets.HEATMAP)
        coverage = get(Datasets.COVERAGE)
        terminations = get(Datasets.TERMINATIONS)
        mi = get(Datasets.MI)
        covariance = get(Datasets.COV)
        norm = get(Datasets.NORM)
        pairwise_snr = get(Datasets.PAIRWISE_SNR)

        # Verify required fields exist
        if reactivity is None:
            raise ValueError(f"Required dataset '{Datasets.REACTIVITY}' not found in file")
        if reads is None:
            raise ValueError(f"Required dataset '{Datasets.READS}' not found in file")
        if error is None:
            raise ValueError(f"Required dataset '{Datasets.ERROR}' not found in file")
        if snr is None:
            raise ValueError(f"Required dataset '{Datasets.SNR}' not found in file")

        n = reactivity.shape[0]
        length = reactivity.shape[1] if reactivity.ndim > 1 else 0

        # Create ProbingData with defaults for missing optional fields
        data = cls(
            sequences=sequences,
            reactivity=reactivity,
            reads=reads,
            error=error,
            snr=snr,
            mask=np.ones_like(reactivity, dtype=bool),
            heatmap=heatmap if heatmap is not None else np.zeros((4, 7)),
            coverage=coverage if coverage is not None else (~np.isnan(reactivity)).astype(float),
            terminations=terminations if terminations is not None else np.zeros((n, length)),
            pairs=None,
            probability=None,
            covariance=covariance,
            mi=mi,
            pairwise_snr=pairwise_snr,
        )
        data.norm = norm
        return data


def _sequences_from_counts(
    file: h5py.File,
    dataset: str = "sequences",
) -> Union[da.Array, None]:
    # Get the tokenized sequences if they exist

    if dataset in file:
        return da.from_array(file[dataset], chunks="auto")

    return None


def _probing_data_from_counts(
    file: h5py.File, groups: DataGroups, opts: Opts, lengths: da.Array
) -> Union[ProbingData, None]:
    if not groups.oneD:
        return None

    sequences = _sequences_from_counts(file)

    counts = _accumulate_arrays(file, groups.oneD)
    if _has_arrays(file, groups.twoD):
        pairs = _accumulate_arrays(file, groups.twoD)
    else:
        pairs = None

    data = _data_from_counts(counts, pairs, opts, lengths)
    # Type ignore: _data_from_counts returns non-None for required fields
    return ProbingData(sequences, *data)  # type: ignore[arg-type]


def _merge_conditions(
    mod: ProbingData,
    nomod: Union[ProbingData, None],
) -> ProbingData:
    if nomod is None:
        return mod

    sequences = mod.sequences

    # 1D data

    reactivity = mod.reactivity - nomod.reactivity
    reads = mod.reads + nomod.reads
    error = np.sqrt(mod.error**2 + nomod.error**2)
    snr = _snr(reactivity, error)
    mask = mod.mask & nomod.mask
    heatmap = mod.heatmap
    coverage = mod.coverage + nomod.coverage
    terminations = mod.terminations + nomod.terminations

    # 2D data

    pairs = mod.pairs
    probability = mod.probability
    mi = mod.mi
    covariance = mod.covariance
    pairwise_snr = mod.pairwise_snr

    return ProbingData(
        sequences,
        reactivity,
        reads,
        error,
        snr,
        mask,
        heatmap,
        coverage,
        terminations,
        pairs,
        probability,
        covariance,
        mi,
        pairwise_snr,
    )


def _lengths_from_fasta(file: str) -> da.Array:
    lengths = []
    length = 0
    start = True

    with open(file) as f:
        for line in f:
            if line.startswith(">"):
                if start:
                    start = False
                    continue
                lengths.append(length)
                length = 0
            else:
                length += len(line.strip())

    lengths.append(length)
    return da.array(lengths)


# Canonical base -> token map for sequence-aware normalization (T folds onto U).
_BASE_TOKEN = {"A": 0, "C": 1, "G": 2, "U": 3, "T": 3}


def _sequences_from_fasta(fasta: str, width: int) -> np.ndarray:
    """Per-position base tokens for each reference (FASTA order).

    A/C/G/U(/T) map to 0/1/2/3; any other character and padding map to -1. Each
    row is padded or truncated to ``width`` so it aligns with the reactivity
    arrays. Used to give sequence-aware norm schemes (e.g. per-nucleotide DMS)
    the base identity when the counts file carries no tokenized sequence.
    """
    seqs: list[str] = []
    cur: list[str] = []
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))

    out = np.full((len(seqs), width), -1, dtype=np.int8)
    for i, s in enumerate(seqs):
        for j, ch in enumerate(s[:width]):
            out[i, j] = _BASE_TOKEN.get(ch.upper(), -1)
    return out


def _names_from_fasta(fasta: str) -> list[str]:
    """Reference names in FASTA order: the first whitespace-delimited token of each
    header (matching the reference names used by the aligner / BAM)."""
    names: list[str] = []
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip()
                names.append(header.split()[0] if header else "")
    return names


def _ensure_sequences(data: ProbingData, fasta: str, opts: Opts) -> None:
    """Attach per-position base tokens from the FASTA when the counts file provided
    none (no-``--tokenize`` case).

    Always attached (not gated on the norm scheme) so the file-level ``sequence``
    dataset is present for every run -- sequence-aware norm schemes and the report
    (which colors each base by its reactivity) both read it. ``opts`` is retained
    for signature stability."""
    del opts  # no longer used; sequences are attached unconditionally
    if data.sequences is None:
        width = int(np.asarray(data.reactivity).shape[1])
        data.sequences = _sequences_from_fasta(fasta, width)


def compute_reactivity(
    file: h5py.File,
    fasta: str,
    opts: Opts,
) -> tuple[ProbingData, Union[ProbingData, None], ProbingData]:
    """
    Compute reactivity profiles from mutation count data.

    Args:
        file: Open HDF5 file containing count data
        fasta: Path to FASTA file with reference sequences
        opts: Processing options

    Returns:
        Tuple of (modified, unmodified, combined) ProbingData objects

    Raises:
        FileNotFoundError: If FASTA file does not exist
    """
    if not os.path.exists(fasta):
        raise FileNotFoundError(f"FASTA file not found: {fasta}")

    # Get reactivity and associated values

    lengths = _lengths_from_fasta(fasta)
    mod = _probing_data_from_counts(file, opts.mod, opts, lengths)
    nomod = _probing_data_from_counts(file, opts.nomod, opts, lengths)

    assert mod is not None
    combined = _merge_conditions(mod, nomod)

    # Compute

    mod = mod.compute()
    combined = combined.compute()
    if nomod is not None:
        nomod = nomod.compute()

    # Normalize

    _ensure_sequences(combined, fasta, opts)
    norm = normalization([combined], opts)[0]

    mod.normalize(norm)
    combined.normalize(norm)
    if nomod is not None:
        nomod.normalize(norm)

    # Clip

    mod.clip(opts.clip[0], opts.clip[1])
    combined.clip(opts.clip[0], opts.clip[1])
    if nomod is not None:
        nomod.clip(opts.clip[0], opts.clip[1])

    return mod, nomod, combined


@dataclass
class ExperimentResult:
    """Per-experiment output of :func:`compute_reactivities`."""

    experiment: Experiment
    mod: ProbingData
    nomod: Union[ProbingData, None]
    combined: ProbingData


def compute_reactivities(
    file: h5py.File,
    fasta: str,
    experiments: list[Experiment],
    opts: Opts,
) -> list[ExperimentResult]:
    """Compute reactivity profiles for one or more experiments.

    All experiments share a single FASTA reference. Mutation counts are read
    once, then per-experiment raw reactivities are computed and normalized. The
    normalization granularity comes from ``opts`` (``per_experiment``,
    ``per_reference``): by default one factor is pooled across all experiments
    and references, keeping the experiments on a comparable scale.

    Args:
        file: Open HDF5 file containing count data for all experiments.
        fasta: Path to FASTA file with reference sequences.
        experiments: One or more named experiments, each with mod and optional
            nomod paths.
        opts: Processing options. ``opts.mod`` and ``opts.nomod`` are ignored;
            data locations come from ``experiments``. The other fields (cutoff,
            ins, dels, norm, blank, clip, sig, and the granularity axes) are
            applied to every experiment.

    Returns:
        One :class:`ExperimentResult` per input experiment, in the same order.

    Raises:
        FileNotFoundError: If the FASTA file does not exist.
        ValueError: If ``experiments`` is empty.
    """
    if not experiments:
        raise ValueError("compute_reactivities requires at least one Experiment.")
    if not os.path.exists(fasta):
        raise FileNotFoundError(f"FASTA file not found: {fasta}")

    lengths = _lengths_from_fasta(fasta)

    # First pass: compute raw (unnormalized, unclipped) reactivities per experiment.

    raw: list[ExperimentResult] = []
    for experiment in experiments:
        mod_groups = DataGroups(list(experiment.mod))
        nomod_groups = DataGroups(list(experiment.nomod) if experiment.nomod else None)

        mod = _probing_data_from_counts(file, mod_groups, opts, lengths)
        nomod = _probing_data_from_counts(file, nomod_groups, opts, lengths)

        if mod is None:
            raise ValueError(
                f"Experiment '{experiment.name}' has no mod datasets; at least one is required."
            )

        combined = _merge_conditions(mod, nomod)

        mod = mod.compute()
        combined = combined.compute()
        if nomod is not None:
            nomod = nomod.compute()

        _ensure_sequences(combined, fasta, opts)
        raw.append(ExperimentResult(experiment=experiment, mod=mod, nomod=nomod, combined=combined))

    # Compute the per-experiment normalization factors (granularity from opts),
    # then apply and clip.

    factors = normalization([r.combined for r in raw], opts)

    for r, norm in zip(raw, factors):
        r.mod.normalize(norm)
        r.combined.normalize(norm)
        if r.nomod is not None:
            r.nomod.normalize(norm)

        r.mod.clip(opts.clip[0], opts.clip[1])
        r.combined.clip(opts.clip[0], opts.clip[1])
        if r.nomod is not None:
            r.nomod.clip(opts.clip[0], opts.clip[1])

    return raw


# --- SNR-vs-read-depth projection ------------------------------------------
# Parameters of the relative-depth axis and the pareto frontier search.
_SNR_DEPTH_MIN = 0.1  # relative-depth axis lower bound (0.1x current depth)
_SNR_DEPTH_MAX = 10.0  # relative-depth axis upper bound (10x current depth)
_SNR_DEPTH_POINTS = 1000  # samples along the relative-depth axis
_SNR_PARETO_SPLITS = 200  # mod/nomod read-allocation fractions for the pareto frontier
_SNR_DEPTH_CHUNK = 100  # depth points per chunk (bounds the mod/nomod curve memory)
_SNR_PRIOR = 0.001  # Beta prior on the mutation rate; sets the per-position variance floor

_SNR_CURVE_KEYS = ("xi", "mod", "mod_sem", "nomod", "nomod_sem", "pareto", "pareto_sem")


@dataclass
class SNRCurves:
    """Expected mean SNR as a function of relative total read depth.

    Computed once in the data layer (memory-bounded) and saved to the reactivity
    HDF5, so the plot backends render rather than recompute. ``xi`` is the
    relative-depth axis; the rest are mean-SNR curves over it. The ``nomod`` and
    ``pareto`` curves are present only when an unmodified control exists.
    """

    xi: np.ndarray
    mod: np.ndarray
    mod_sem: np.ndarray
    nomod: Union[np.ndarray, None] = None
    nomod_sem: Union[np.ndarray, None] = None
    pareto: Union[np.ndarray, None] = None
    pareto_sem: Union[np.ndarray, None] = None

    def save(self, group: str, path: str) -> None:
        """Append the curves to an existing HDF5 file under ``group/snr-<name>``."""
        with h5py.File(path, "a") as f:
            grp = f if group == "" else f.require_group(group)
            for key in _SNR_CURVE_KEYS:
                value = getattr(self, key)
                if value is not None:
                    grp.create_dataset(f"snr-{key}", data=np.asarray(value))

    @classmethod
    def load(cls, group: str, file: h5py.File) -> Union[SNRCurves, None]:
        prefix = f"{group}/snr-" if group else "snr-"
        if f"{prefix}xi" not in file:
            return None
        vals = {
            k: (np.array(file[f"{prefix}{k}"]) if f"{prefix}{k}" in file else None)
            for k in _SNR_CURVE_KEYS
        }
        return cls(**vals)


def compute_snr_curves(
    mod: ProbingData,
    nomod: Union[ProbingData, None],
    combined: ProbingData,
) -> SNRCurves:
    """Expected mean SNR as a function of relative total read depth.

    Reduces over references/positions at each depth, so the result is small
    ``(_SNR_DEPTH_POINTS,)`` curves (see :class:`SNRCurves`). The pareto search
    is chunked to bound memory.
    """
    reactivity = np.asarray(combined.reactivity)
    mod_err = np.asarray(mod.error)
    mod_reads = float(np.asarray(mod.reads).max())

    prior = _SNR_PRIOR
    mod_err2 = np.maximum(mod_err**2, prior * (1 - prior) / mod_reads)

    if nomod is not None:
        nomod_err = np.asarray(nomod.error)
        nomod_reads = float(np.asarray(nomod.reads).max())
        total_reads = mod_reads + nomod_reads
        nomod_err2: Union[np.ndarray, None] = np.maximum(
            nomod_err**2, prior * (1 - prior) / nomod_reads
        )
    else:
        nomod_err2 = None
        total_reads = mod_reads

    # SEM of mean SNR at current depth; bands at projected depths scale this.
    current_se2 = mod_err2.copy()
    if nomod_err2 is not None:
        current_se2 = current_se2 + nomod_err2
    current_snr = np.where(np.sqrt(current_se2) > 0, reactivity / np.sqrt(current_se2), 0.0)
    valid = current_snr.ravel()[~np.isnan(current_snr.ravel())]
    snr_at_1 = valid.mean() if len(valid) > 0 else 1.0
    snr_sem_at_1 = valid.std() / np.sqrt(len(valid)) if len(valid) > 1 else 0.0

    def mean_snr(mod_scales: np.ndarray, nomod_scales: np.ndarray) -> np.ndarray:
        # Chunk the depth axis so the intermediate ``(depth, refs, positions)``
        # array stays bounded on reference-heavy libraries.
        out = np.empty(len(mod_scales))
        for i in range(0, len(mod_scales), _SNR_DEPTH_CHUNK):
            sl = slice(i, i + _SNR_DEPTH_CHUNK)
            se2 = mod_err2[None] / mod_scales[sl, None, None]
            if nomod_err2 is not None:
                se2 = se2 + nomod_err2[None] / nomod_scales[sl, None, None]
            se = np.sqrt(se2)
            snr = np.where(se > 0, reactivity[None] / se, 0.0)
            out[sl] = np.nanmean(snr, axis=-1).mean(axis=-1)
        return out

    def band(curve: np.ndarray) -> np.ndarray:
        ratio = np.where(snr_at_1 > 0, curve / snr_at_1, 0.0)
        return np.abs(ratio) * snr_sem_at_1

    xi = np.geomspace(_SNR_DEPTH_MIN, _SNR_DEPTH_MAX, _SNR_DEPTH_POINTS)

    if nomod is None:
        snr_mod = mean_snr(xi, np.ones_like(xi))
        return SNRCurves(xi=xi, mod=snr_mod, mod_sem=band(snr_mod))

    # Modified: extra reads go to mod; Unmodified: extra reads go to nomod.
    mod_scales = np.maximum((xi * total_reads - nomod_reads) / mod_reads, 1e-10)
    snr_mod = mean_snr(mod_scales, np.ones_like(mod_scales))
    nomod_scales = np.maximum((xi * total_reads - mod_reads) / nomod_reads, 1e-10)
    snr_nomod = mean_snr(np.ones_like(nomod_scales), nomod_scales)

    # Pareto frontier: best mod/nomod read split at each depth. At a split
    # fraction ``f`` and relative depth ``xi``, ``mean_snr`` factorizes as
    # ``sqrt(xi * total_reads) * G(f)`` with
    # ``G(f) = mean_r nanmean_p [ reactivity / sqrt(a/f + b/(1-f)) ]``,
    # where ``a = mod_err2 * mod_reads`` and ``b = nomod_err2 * nomod_reads``.
    # Depth only rescales the objective by a positive constant, so the optimal
    # split is depth-independent: evaluate ``G`` over the fraction grid once,
    # then scale its maximum across all depths.
    assert nomod_err2 is not None  # guaranteed by the nomod-is-None early return
    fracs = np.linspace(0.01, 0.99, _SNR_PARETO_SPLITS)
    a = mod_err2 * mod_reads
    b = nomod_err2 * nomod_reads
    g_of_f = np.empty(len(fracs))
    for j, f in enumerate(fracs):
        se = np.sqrt(a / f + b / (1 - f))
        snr = np.where(se > 0, reactivity / se, 0.0)
        g_of_f[j] = np.nanmean(snr, axis=-1).mean()
    snr_pareto = np.sqrt(xi * total_reads) * g_of_f.max()

    return SNRCurves(
        xi=xi,
        mod=snr_mod,
        mod_sem=band(snr_mod),
        nomod=snr_nomod,
        nomod_sem=band(snr_nomod),
        pareto=snr_pareto,
        pareto_sem=band(snr_pareto),
    )


def save_groups(
    path: str,
    results: list[tuple[str, ProbingData]],
    names: Union[list[str], None] = None,
) -> None:
    """Save multiple ProbingData objects to a single HDF5 file.

    Each experiment is written under its own group at the file root. Shared,
    experiment-independent datasets (the reference sequences and names) go under a
    single ``meta`` group, leaving room for future metadata.

    Args:
        path: Output HDF5 file. Overwritten if it exists.
        results: List of ``(group_name, ProbingData)`` tuples. ``meta`` is reserved
            for the metadata group and may not be used as an experiment name.
        names: Reference names in FASTA order (see ``_names_from_fasta``). Written
            as the ``meta/name`` dataset when provided, so downstream consumers
            (e.g. the report) can label references by name rather than index.
    """
    if any(name == Datasets.META for name, _ in results):
        raise ValueError(
            f"'{Datasets.META}' is reserved for the metadata group and "
            "cannot be used as an experiment name."
        )
    with h5py.File(path, "w") as f:
        for name, data in results:
            grp = f if name == "" else f.create_group(name)
            for key, arr in data.datasets().items():
                grp.create_dataset(key, data=arr)

        sequences = results[0][1].sequences if results else None
        if sequences is not None or names:
            meta = f.create_group(Datasets.META)
            if sequences is not None:
                meta.create_dataset(Datasets.SEQUENCE, data=np.asarray(sequences))
            if names:
                meta.create_dataset(
                    Datasets.NAME, data=np.array(names, dtype=h5py.string_dtype(encoding="utf-8"))
                )
