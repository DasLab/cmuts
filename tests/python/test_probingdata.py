"""Round-trip tests for ProbingData HDF5 save/load.

The reactivity h5 must be a complete record of the plottable data so that
`cmuts plot` can regenerate every figure `cmuts normalize` makes without the
original counts. In particular coverage and terminations (previously faked on
load) must survive a save/load cycle.
"""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from cmuts.internal import ProbingData, save_groups

pytestmark = pytest.mark.no_external_dependencies


def _make(n: int = 5, length: int = 30) -> ProbingData:
    rng = np.random.default_rng(0)
    data = ProbingData(
        sequences=None,
        reactivity=rng.random((n, length)),
        reads=rng.integers(1, 1000, n).astype(float),
        error=rng.random((n, length)),
        snr=rng.random((n, length)),
        mask=np.ones((n, length), dtype=bool),
        heatmap=rng.random((4, 7)),
        coverage=rng.random(length),
        terminations=rng.random((n, length)),
        pairs=None,
        probability=None,
        covariance=None,
        mi=None,
    )
    data.norm = np.ones(1)
    return data


def test_save_load_preserves_plot_inputs(tmp_path):
    pd = _make()
    out = tmp_path / "reactivity.h5"
    pd.save("grp", str(out))

    with h5py.File(out, "r") as f:
        assert {"coverage", "terminations"} <= set(f["grp"].keys())
        loaded = ProbingData.load("grp", f)

    for field in ("reactivity", "reads", "error", "snr", "heatmap", "coverage", "terminations"):
        assert np.allclose(getattr(loaded, field), getattr(pd, field)), field


def test_save_groups_is_lossless(tmp_path):
    """The multi-group writer used by `cmuts normalize` must also persist
    coverage/terminations so `cmuts plot` can read them back."""
    pd = _make()
    out = tmp_path / "profiles.h5"
    save_groups(str(out), [("g1", pd)])
    with h5py.File(out, "r") as f:
        assert {"coverage", "terminations"} <= set(f["g1"].keys())
        loaded = ProbingData.load("g1", f)
    assert np.allclose(loaded.coverage, pd.coverage)
    assert np.allclose(loaded.terminations, pd.terminations)


def test_load_falls_back_for_legacy_files(tmp_path):
    """An h5 without coverage/terminations (older normalize) still loads."""
    pd = _make()
    out = tmp_path / "legacy.h5"
    pd.save("grp", str(out))
    with h5py.File(out, "a") as f:
        del f["grp/coverage"]
        del f["grp/terminations"]
    with h5py.File(out, "r") as f:
        loaded = ProbingData.load("grp", f)
    # Reconstructed (per reference), not crashed.
    assert loaded.coverage.shape == pd.reactivity.shape
    assert loaded.terminations.shape == pd.reactivity.shape


def test_per_reference_coverage_round_trips(tmp_path):
    """Coverage is stored per reference (n_ref, L) so the report can plot a single
    sequence's coverage; the shape and values must survive save/load."""
    pd = _make(n=4, length=20)
    rng = np.random.default_rng(1)
    pd.coverage = rng.random((4, 20))  # per-reference, as compute_reactivities now writes
    out = tmp_path / "profiles.h5"
    save_groups(str(out), [("g1", pd)])
    with h5py.File(out, "r") as f:
        loaded = ProbingData.load("g1", f)
    assert loaded.coverage.shape == (4, 20)
    assert np.allclose(loaded.coverage, pd.coverage)


def test_save_groups_writes_reference_names(tmp_path):
    """Reference names (FASTA order) are persisted under meta/name so the report
    can label sequences by name rather than index."""
    pd = _make()
    out = tmp_path / "profiles.h5"
    save_groups(str(out), [("g1", pd)], names=["tRNA", "16S", "riboswitch"])
    with h5py.File(out, "r") as f:
        assert "name" in f["meta"]
        names = [n.decode() if isinstance(n, bytes) else n for n in f["meta"]["name"]]
    assert names == ["tRNA", "16S", "riboswitch"]


def test_save_groups_writes_sequence_under_meta(tmp_path):
    """Sequences (base tokens) are persisted under meta/sequence and round-trip."""
    pd = _make(n=3, length=10)
    pd.sequences = np.tile(np.arange(10, dtype=np.int8), (3, 1))
    out = tmp_path / "profiles.h5"
    save_groups(str(out), [("g1", pd)])
    with h5py.File(out, "r") as f:
        assert "sequence" in f["meta"] and "sequence" not in f
        loaded = ProbingData.load("g1", f)
    assert np.array_equal(loaded.sequences, pd.sequences)


def test_save_groups_rejects_meta_experiment_name(tmp_path):
    """'meta' is reserved for the metadata group, not an experiment name."""
    pd = _make()
    with pytest.raises(ValueError, match="reserved"):
        save_groups(str(tmp_path / "x.h5"), [("meta", pd)])
