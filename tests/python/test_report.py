"""Contract tests for the HDF5-driven HTML report (``cmuts.report.build``).

These assert *behavior* -- which plots and data the report includes, and when --
rather than presentation (headings, CSS classes, colors, layout), which changes
freely. The stable anchors are the client-JS contract identifiers: the embedded
``cmuts-figs`` / ``cmuts-data`` JSON blobs, the ``seqselect`` dropdown, and the
``multiplot`` / ``diffplot`` targets.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import h5py
import numpy as np
import pytest

pytest.importorskip("plotly.graph_objects")

from cmuts.report import build  # noqa: E402

pytestmark = pytest.mark.no_external_dependencies


def _write_h5(
    path: Path,
    groups: list,
    n_ref: int = 3,
    length: int = 12,
    *,
    matrices: bool = True,
    snr: bool = True,
    coverage_2d: bool = True,
    names: bool = True,
) -> None:
    """Write a synthetic reactivity h5 using the normalize<->plot dataset names."""
    row = np.linspace(0.05, 0.95, length, dtype="f4")
    with h5py.File(path, "w") as f:
        for g in groups:
            grp = f.create_group(g)
            grp["reactivity"] = np.tile(row, (n_ref, 1))
            grp["error"] = np.full((n_ref, length), 0.05, dtype="f4")
            grp["reads"] = (np.arange(1, n_ref + 1) * 100).astype("f4")
            grp["SNR"] = np.full(n_ref, 2.0, dtype="f4")
            grp["heatmap"] = np.full((4, 7), 0.1, dtype="f4")
            grp["coverage"] = (
                np.tile(row, (n_ref, 1)) * 100 if coverage_2d else row * 100
            ).astype("f4")
            grp["terminations"] = np.abs(np.tile(row, (n_ref, 1))).astype("f4")
            if matrices:
                base = np.eye(length, dtype="f4")[None]
                grp["mutual-information"] = np.repeat(base * 0.1, n_ref, axis=0)
                grp["covariance"] = np.repeat(base, n_ref, axis=0)
            if snr:
                grp["snr-xi"] = np.geomspace(0.1, 10, 1000)
                grp["snr-mod"] = np.linspace(0, 2, 1000)
                grp["snr-mod_sem"] = np.full(1000, 0.1)
        # cmuts normalize always writes meta/sequence; meta/name only with a FASTA.
        meta = f.create_group("meta")
        meta.create_dataset(  # ACGU tokens (0-3) cycled across each reference
            "sequence", data=np.tile(np.arange(length, dtype=np.int8) % 4, (n_ref, 1)))
        if names:
            meta.create_dataset(
                "name",
                data=np.array([f"ref{i}" for i in range(n_ref)],
                              dtype=h5py.string_dtype(encoding="utf-8")),
            )


def _embed(html: str, tag_id: str):
    """Parse an embedded ``application/json`` blob by id, or None if absent."""
    m = re.search(
        r'<script type="application/json" id="' + tag_id + r'">(.*?)</script>',
        html, re.DOTALL,
    )
    return json.loads(m.group(1)) if m else None


# --- default vs --all embedding --------------------------------------------

def test_default_mode_embeds_no_per_reference_data(tmp_path: Path) -> None:
    html = build(str(_h5(tmp_path, ["A", "B"], n_ref=3)))
    assert _embed(html, "cmuts-figs") is None
    assert _embed(html, "cmuts-data") is None
    assert 'class="seqselect"' not in html      # no per-sequence dropdown
    assert 'id="multiplot"' not in html         # no comparison sections
    assert 'id="diffplot"' not in html
    assert "plotly-graph-div" in html           # aggregate plots still rendered


def test_all_mode_embeds_one_entry_per_reference(tmp_path: Path) -> None:
    html = build(str(_h5(tmp_path, ["A", "B"], n_ref=3)), embed_all=True)
    figs = _embed(html, "cmuts-figs")
    assert set(figs) == {"A", "B"}
    for group in figs.values():
        assert set(group) == {"0", "1", "2"}   # one embedded figure set per reference
        # The profile is client-rendered from the raw arrays, not pre-built here.
        assert set(group["0"]) == {"coverage", "terminations", "mi", "correlation"}
    # Profiles are drawn from the embedded arrays instead.
    assert _embed(html, "cmuts-data")["A"]["reactivity"]
    assert 'class="seqselect"' in html


def test_missing_matrices_are_omitted(tmp_path: Path) -> None:
    figs = _embed(
        build(str(_h5(tmp_path, ["A"], n_ref=2, matrices=False)), embed_all=True),
        "cmuts-figs",
    )
    assert set(figs["A"]["0"]) == {"coverage", "terminations"}


def test_legacy_1d_coverage_omits_per_sequence_coverage(tmp_path: Path) -> None:
    # Aggregate coverage still works from a 1-D array, but a single reference's
    # coverage can only be shown when coverage is stored per reference (2-D).
    figs = _embed(
        build(str(_h5(tmp_path, ["A"], n_ref=2, coverage_2d=False)), embed_all=True),
        "cmuts-figs",
    )
    assert "coverage" not in figs["A"]["0"]


# --- single-reference layout ------------------------------------------------

def test_single_reference_embeds_profile_without_all(tmp_path: Path) -> None:
    html = build(str(_h5(tmp_path, ["A"], n_ref=1)))   # note: no embed_all
    assert 'class="seqselect"' not in html             # nothing to select
    assert "js-profile" in html                        # profile is client-rendered
    assert _embed(html, "cmuts-data")["A"]["reactivity"]  # its arrays are embedded
    assert "plotly-graph-div" in html                  # the other plots render inline


# --- comparison-section gating (multi-plot / difference) --------------------

def test_no_comparison_for_single_experiment_single_reference(tmp_path: Path) -> None:
    html = build(str(_h5(tmp_path, ["A"], n_ref=1)), embed_all=True)
    assert 'id="multiplot"' not in html and 'id="diffplot"' not in html


def test_comparison_for_single_experiment_multiple_references(tmp_path: Path) -> None:
    html = build(str(_h5(tmp_path, ["A"], n_ref=3)), embed_all=True)
    assert 'id="multiplot"' in html and 'id="diffplot"' in html


def test_comparison_for_multiple_experiments(tmp_path: Path) -> None:
    html = build(str(_h5(tmp_path, ["A", "B"], n_ref=1)), embed_all=True)
    assert 'id="multiplot"' in html and 'id="diffplot"' in html


# --- embedded data for the multi-plot / difference JS ----------------------

def test_data_embed_has_arrays_labels_and_stats(tmp_path: Path) -> None:
    data = _embed(build(str(_h5(tmp_path, ["A"], n_ref=2)), embed_all=True), "cmuts-data")
    a = data["A"]
    assert a["n"] == 2
    assert len(a["reactivity"]) == 2 and len(a["reactivity"][0]) == 12
    assert len(a["error"]) == 2
    assert a["names"] == ["ref0", "ref1"]              # FASTA names used
    assert len(a["seqstats"]) == 2 and "reads" in a["seqstats"][0]


def test_reference_labels_fall_back_without_names(tmp_path: Path) -> None:
    data = _embed(
        build(str(_h5(tmp_path, ["A"], n_ref=2, names=False)), embed_all=True), "cmuts-data"
    )
    assert data["A"]["names"] == ["Reference 1", "Reference 2"]


def test_embedded_arrays_are_json_safe(tmp_path: Path) -> None:
    # NaNs must be emitted as JSON null, or the browser's JSON.parse would fail.
    path = tmp_path / "nan.h5"
    with h5py.File(path, "w") as f:
        g = f.create_group("A")
        react = np.full((1, 5), 0.5, dtype="f4")
        react[0, 2] = np.nan
        g["reactivity"] = react
        g["error"] = np.full((1, 5), 0.05, dtype="f4")
        g["reads"] = np.array([100.0], dtype="f4")
        g["SNR"] = np.array([2.0], dtype="f4")
        g["heatmap"] = np.full((4, 7), 0.1, dtype="f4")
        g["coverage"] = np.full((1, 5), 100.0, dtype="f4")
        g["terminations"] = np.full((1, 5), 0.1, dtype="f4")
    html = build(str(path), embed_all=True)
    blob = re.search(r'id="cmuts-data">(.*?)</script>', html, re.DOTALL).group(1)
    assert "NaN" not in blob and "null" in blob


# --- filtering & structure --------------------------------------------------

def test_group_filter_restricts_to_one_group(tmp_path: Path) -> None:
    data = _embed(
        build(str(_h5(tmp_path, ["A", "B"], n_ref=1)), embed_all=True, group="A"), "cmuts-data"
    )
    assert set(data) == {"A"}


def test_structure_block_is_included(tmp_path: Path) -> None:
    html = build(str(_h5(tmp_path, ["A"], n_ref=1)), structure="color /A red")
    assert "color /A red" in html


# --- colored sequence view --------------------------------------------------

# The colored sequence view is client-rendered from the cmuts-seq embed blob (the
# per-reference base strings), so these assert on that contract, not on markup.

def test_sequence_embedded_per_reference(tmp_path: Path) -> None:
    seq = _embed(build(str(_h5(tmp_path, ["A"], n_ref=3)), embed_all=True), "cmuts-seq")
    assert seq is not None
    assert len(seq) == 3 and all(len(s) == 12 for s in seq)   # one base string per reference


def test_sequence_tokens_decode_to_bases(tmp_path: Path) -> None:
    # Even without --all, a single reference embeds its arrays and sequence; tokens
    # 0/1/2/3 (cycled) must decode to ACGU.
    seq = _embed(build(str(_h5(tmp_path, ["A"], n_ref=1))), "cmuts-seq")
    assert seq == ["ACGUACGUACGU"]


def _h5(tmp_path: Path, groups: list, **kw) -> Path:
    path = tmp_path / "profiles.h5"
    _write_h5(path, groups, **kw)
    return path
