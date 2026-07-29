"""Plotly implementations of cmuts plots.

Each function returns a ``go.Figure``. Every builder takes the raw arrays it
needs (read from the reactivity HDF5 by the caller), so this module depends only
on numpy/plotly and never on the normalize pipeline internals.
"""

from typing import Any, Optional

import numpy as np
import plotly.colors as pc
import plotly.graph_objects as go

from . import _transforms

# ===========================================================================
# Styling constants
# ===========================================================================
FONT_FAMILY = "Helvetica"
LABEL_SIZE = 14
TITLE_SIZE = 15
TICK_SIZE = 13
LEGEND_SIZE = 12

# Tight margins so the data area fills its container (the report figure card).
_MARGIN = {"l": 55, "r": 20, "t": 40, "b": 45}

# Pre-sampled colours: a light fill + darker line for area plots, and a mid tone
# for the SNR confidence bands (modified = RdPu, unmodified = PuBu).
_RDPU_FILL = pc.sample_colorscale("RdPu", [0.3])[0]
_RDPU_LINE = pc.sample_colorscale("RdPu", [0.8])[0]
_RDPU_MID = pc.sample_colorscale("RdPu", [0.7])[0]
_PUBU_MID = pc.sample_colorscale("PuBu", [0.7])[0]


# ===========================================================================
# Shared trace / layout / colorbar helpers
# ===========================================================================
def _title(base: str, name: str) -> str:
    """`"Base (name)"` when a group name is given, else just `"Base"`."""
    return f"{base} ({name})" if name else base


def _rgba(rgb_str: str, alpha: float) -> str:
    return rgb_str.replace("rgb(", "rgba(").replace(")", f", {alpha})")


def _base_layout(title: str, xlabel: str = "", ylabel: str = "") -> dict[str, Any]:
    return {
        "font": {"family": FONT_FAMILY},
        "title": {"text": title, "font": {"size": TITLE_SIZE}},
        "xaxis": {
            "title": {"text": xlabel, "font": {"size": LABEL_SIZE}},
            "tickfont": {"size": TICK_SIZE},
        },
        "yaxis": {
            "title": {"text": ylabel, "font": {"size": LABEL_SIZE}},
            "tickfont": {"size": TICK_SIZE},
            "gridcolor": "rgba(0,0,0,0.15)",
        },
        "plot_bgcolor": "white",
        "legend": {"font": {"size": LEGEND_SIZE}},
        "margin": _MARGIN,
    }


def _colorbar(title: str, **extra: Any) -> dict[str, Any]:
    """A colorbar styled like the old matplotlib figures: outlined, with outside
    ticks and a rotated title on the right. ``extra`` overrides/adds keys (e.g.
    tickvals/ticktext). Length is tied to the plotting area (see the heatmap
    builders' xaxis domain) so it matches the plot height."""
    cb = {
        "title": {"text": title, "side": "right", "font": {"size": LABEL_SIZE}},
        "outlinecolor": "black",
        "outlinewidth": 1,
        "ticks": "outside",
        "ticklen": 5,
        "tickcolor": "black",
        "tickfont": {"size": TICK_SIZE},
        "thickness": 15,
    }
    cb.update(extra)
    return cb


def _line_trace(
    fig: go.Figure,
    data: np.ndarray,
    fill_color: str = _RDPU_FILL,
    line_color: str = _RDPU_LINE,
    label: str = "",
) -> None:
    data = np.asarray(data)
    x = np.arange(data.shape[0])
    fig.add_trace(
        go.Scatter(
            x=x,
            y=data,
            fill="tozeroy",
            fillcolor=_rgba(fill_color, 0.5),
            line={"color": line_color, "width": 1},
            name=label,
            showlegend=bool(label),
        )
    )


def _add_band(
    fig: go.Figure,
    x: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    fillcolor: str,
) -> None:
    """Shade the region between ``lower`` and ``upper``: an invisible lower
    trace followed by an upper trace that fills down to it (``fill="tonexty"``).
    """
    fig.add_trace(go.Scatter(x=x, y=lower, line={"width": 0}, showlegend=False, hoverinfo="skip"))
    fig.add_trace(
        go.Scatter(
            x=x,
            y=upper,
            fill="tonexty",
            fillcolor=fillcolor,
            line={"width": 0},
            showlegend=False,
            hoverinfo="skip",
        )
    )


def _line_figure(
    data: np.ndarray,
    title: str,
    xlabel: str,
    ylabel: str,
    *,
    yrange: Optional[list] = None,
    axis_lines: bool = False,
) -> go.Figure:
    """A single fill-to-zero line plot -- the shape shared by coverage,
    terminations, reads-per-block, etc. ``axis_lines`` draws black x/y spines and
    ``yrange`` fixes the y extent."""
    fig = go.Figure()
    _line_trace(fig, data)
    fig.update_layout(**_base_layout(title, xlabel, ylabel))
    spine = {"showline": True, "linecolor": "black", "linewidth": 1} if axis_lines else {}
    fig.update_xaxes(range=[0, len(data)], **spine)
    if yrange is not None:
        fig.update_yaxes(range=yrange, **spine)
    elif spine:
        fig.update_yaxes(**spine)
    return fig


def _square_matrix_layout(fig: go.Figure) -> None:
    """Square, reversed-y layout for an L x L matrix heatmap. scaleanchor keeps
    the cells square; the wide ``.square`` report card leaves the *width* with
    slack, so the x-domain shrinks to square the data while the y-axis keeps the
    full height -- so the colorbar (full height) matches the plot exactly."""
    fig.update_layout(
        yaxis_autorange="reversed",
        yaxis_scaleanchor="x",
        xaxis_constrain="domain",
        yaxis_constrain="domain",
    )


def _matrix_heatmap(
    z: np.ndarray,
    name: str,
    base_title: str,
    colorscale: Any,
    zmin: float,
    zmax: float,
    colorbar: dict,
) -> go.Figure:
    """An L x L matrix as a square, reversed-y heatmap with a styled colorbar."""
    fig = go.Figure(
        data=go.Heatmap(z=z, colorscale=colorscale, zmin=zmin, zmax=zmax, colorbar=colorbar)
    )
    fig.update_layout(**_base_layout(_title(base_title, name), "Residue", "Residue"))
    _square_matrix_layout(fig)
    return fig


# ===========================================================================
# Per-base mutation heatmap
# ===========================================================================
_HEATMAP_NTS = ["A", "C", "G", "U"]
_HEATMAP_MODS = ["A", "C", "G", "U", "del", "ins", "term"]


def plot_heatmap(heatmap: np.ndarray, name: str = "") -> go.Figure:
    heatmap = np.asarray(heatmap)
    with np.errstate(divide="ignore"):
        log_heatmap = np.where(heatmap > 0, np.log10(heatmap), np.nan)

    # Build descriptive hover text per cell
    hover_text = []
    for i, nt in enumerate(_HEATMAP_NTS):
        row = []
        for j, mod in enumerate(_HEATMAP_MODS):
            val = heatmap[i, j]
            prob = f"{val:.4e}" if val > 0 else "0"
            if mod in _HEATMAP_NTS and mod == nt:
                row.append(f"Match ({nt})<br>Probability: {prob}")
            elif mod in _HEATMAP_NTS:
                row.append(f"Mismatch {nt} \u2192 {mod}<br>Probability: {prob}")
            elif mod == "del":
                row.append(f"Deletion of {nt}<br>Probability: {prob}")
            elif mod == "ins":
                row.append(f"Insertion at {nt}<br>Probability: {prob}")
            else:
                row.append(f"Termination at {nt}<br>Probability: {prob}")
        hover_text.append(row)

    title = _title("Modification Heatmap", name)
    fig = go.Figure(
        data=go.Heatmap(
            z=log_heatmap,
            x=_HEATMAP_MODS,
            y=_HEATMAP_NTS,
            colorscale="RdPu",
            zmin=-4,
            zmax=0,
            text=hover_text,
            hoverinfo="text",
            colorbar={
                **_colorbar("Probability"),
                "tickvals": [-4, -3, -2, -1, 0],
                "ticktext": [
                    "10\u207b\u2074",
                    "10\u207b\u00b3",
                    "10\u207b\u00b2",
                    "10\u207b\u00b9",
                    "10\u2070",
                ],
            },
        )
    )

    # Cell outlines
    for i in range(len(_HEATMAP_NTS)):
        for j in range(len(_HEATMAP_MODS)):
            fig.add_shape(
                type="rect",
                x0=j - 0.5,
                x1=j + 0.5,
                y0=i - 0.5,
                y1=i + 0.5,
                line={"color": "black", "width": 1},
                layer="above",
            )

    # Not using _base_layout: this heatmap has a reversed y-axis and its own
    # categorical axis titles. Square cells (scaleanchor + constrain "domain")
    # give the data portion a 7x4 = 7:4 aspect; no fixed width/height, so the
    # figure stays responsive and fits (centers) within its report card. Keep in
    # sync with _base_layout's font sizing by hand if the base style changes.
    fig.update_layout(
        font={"family": FONT_FAMILY},
        title={"text": title, "font": {"size": TITLE_SIZE}},
        xaxis_title="Modification Type",
        xaxis={"showgrid": False, "zeroline": False, "constrain": "domain"},
        yaxis_title="Reference Nucleotide",
        yaxis={
            "autorange": "reversed",
            "showgrid": False,
            "zeroline": False,
            "ticklabelstandoff": 10,
            "scaleanchor": "x",
            "constrain": "domain",
        },
        template="plotly_white",
        margin=_MARGIN,
    )
    return fig


# ===========================================================================
# Position (line) plots
# ===========================================================================
def plot_read_hist(reads: np.ndarray, name: str = "") -> go.Figure:
    bin_centers, counts, normed = _transforms.read_histogram(reads)
    colors = [pc.sample_colorscale("RdPu", [float(v)])[0] for v in normed]

    title = _title("Read distribution", name)
    fig = go.Figure(
        data=go.Bar(
            x=bin_centers,
            y=counts,
            marker_color=colors,
            showlegend=False,
        )
    )
    fig.update_layout(**_base_layout(title, "log10 Read Depth", "Count"))
    fig.update_layout(bargap=0)
    return fig


def plot_reads_per_block(reads: np.ndarray, name: str = "", nblocks: int = 100) -> go.Figure:
    """Mean reads per reference, references split into ``nblocks`` equal-sized
    blocks in FASTA order (no sorting, so spatial trends are preserved)."""
    means = _transforms.reads_per_block(reads, nblocks)
    return _line_figure(
        means, _title("Mean reads per reference bin", name), "Reference bin", "Mean reads"
    )


def plot_termination(term: np.ndarray, name: str = "") -> go.Figure:
    """Termination density for one reference's raw 1-D termination counts.

    Row-normalizes the input to a density. Passing an already-normalized density
    (e.g. the cross-reference aggregate from ``_transforms.termination_density``)
    is a no-op, so callers plot a single reference with ``terminations[i]`` and
    the aggregate with ``termination_density(terminations)``.
    """
    term = np.asarray(term, dtype=float)
    total = term.sum()
    if total > 0:
        term = term / total
    return _line_figure(
        term,
        _title("Termination by position", name),
        "Residue",
        "Termination density",
        yrange=[0, 1.05],
        axis_lines=True,
    )


def plot_coverage(coverage: np.ndarray, reads: np.ndarray, name: str = "") -> go.Figure:
    data = _transforms.coverage_fraction(coverage, reads)
    return _line_figure(
        data,
        _title("Coverage by position", name),
        "Residue",
        "Fraction of reads",
        yrange=[0, 1.05],
        axis_lines=True,
    )


# ===========================================================================
# Matrix heatmaps (reference x position, and L x L pairwise)
# ===========================================================================
def plot_reactivity_heatmap(reactivity: np.ndarray, name: str = "", num: int = 250) -> go.Figure:
    """Reactivity of every reference as a (reference x position) heatmap.

    Clipped to the first ``num`` references so it stays bounded (and renders
    without ``--all``) even for reference-heavy experiments.
    """
    sentinel = -0.01
    data = np.asarray(reactivity)[:num].copy()
    data[np.isnan(data)] = sentinel

    # Grey for the NaN sentinel just below zero, then the RdPu ramp above it.
    eps = 0.001
    colorscale = [
        [0.0, "grey"],
        [eps, "grey"],
        [eps + 0.001, pc.sample_colorscale("RdPu", [0.0])[0]],
        [1.0, pc.sample_colorscale("RdPu", [1.0])[0]],
    ]

    fig = go.Figure(
        data=go.Heatmap(
            z=data,
            colorscale=colorscale,
            zmin=sentinel,
            zmax=1,
            colorbar=_colorbar("Reactivity"),
        )
    )
    fig.update_layout(
        **_base_layout(_title("Heatmap of profiles", name), "Residue", "Sequence Index")
    )
    fig.update_layout(yaxis_autorange="reversed")
    return fig


def plot_pairwise_coverage(values: np.ndarray, name: str = "") -> go.Figure:
    vlow = _transforms.pairwise_coverage_bounds(values).vlow
    with np.errstate(divide="ignore"):
        z = np.log10(np.clip(values, vlow, 1.0))
    return _matrix_heatmap(
        z, name, "Pairwise Coverage", "RdPu", np.log10(vlow), 0, _colorbar("Pairwise Coverage")
    )


def plot_correlation(values: np.ndarray, name: str = "") -> go.Figure:
    linthresh = _transforms.CORRELATION_LINTHRESH
    z = _transforms.symlog(values, linthresh)
    t_max = float(_transforms.symlog(np.array([1.0]), linthresh)[0])
    colorbar = {
        **_colorbar("Pairwise Correlation"),
        "tickvals": [
            _transforms.symlog(np.array([v]), linthresh)[0]
            for v in [-1, -0.1, -0.01, 0, 0.01, 0.1, 1]
        ],
        "ticktext": ["-1", "-0.1", "-0.01", "0", "0.01", "0.1", "1"],
    }
    return _matrix_heatmap(z, name, "Correlation", "PiYG", -t_max, t_max, colorbar)


def plot_mi(values: np.ndarray, name: str = "") -> go.Figure:
    vlow, vhigh = _transforms.mi_bounds(values)
    with np.errstate(divide="ignore"):
        z = np.log10(np.clip(values, vlow, vhigh))
    return _matrix_heatmap(
        z,
        name,
        "Mutual Information",
        "RdPu",
        np.log10(vlow),
        np.log10(vhigh),
        _colorbar("Mutual information"),
    )


# ===========================================================================
# SNR vs read depth
# ===========================================================================
def plot_snr_scaling(
    xi: np.ndarray,
    mod: np.ndarray,
    mod_sem: np.ndarray,
    nomod: Optional[np.ndarray] = None,
    nomod_sem: Optional[np.ndarray] = None,
    pareto: Optional[np.ndarray] = None,
    pareto_sem: Optional[np.ndarray] = None,
    name: str = "",
) -> go.Figure:
    """Render precomputed SNR-vs-read-depth curves (see cmuts.compute_snr_curves).

    ``xi`` is the relative-depth axis; ``mod``/``mod_sem`` the modified mean-SNR
    curve and its band. The ``nomod`` and ``pareto`` curves (with their SEMs) are
    present only when an unmodified control exists; they are drawn together.
    """

    def _trim(snr: np.ndarray) -> slice:
        nz = np.nonzero(snr > 0.01)[0]
        start = max(nz[0] - 1, 0) if len(nz) > 0 else 0
        return slice(start, None)

    fig = go.Figure()

    if nomod is not None:
        # nomod_sem, pareto and pareto_sem are populated together with nomod.
        assert nomod_sem is not None
        assert pareto is not None and pareto_sem is not None
        s = _trim(mod)
        _add_band(
            fig,
            xi[s],
            (mod - mod_sem)[s],
            (mod + mod_sem)[s],
            _rgba(_RDPU_MID, 0.2),
        )
        fig.add_trace(
            go.Scatter(x=xi[s], y=mod[s], line={"color": _RDPU_MID, "width": 2}, name="Modified")
        )

        s = _trim(nomod)
        _add_band(
            fig,
            xi[s],
            (nomod - nomod_sem)[s],
            (nomod + nomod_sem)[s],
            _rgba(_PUBU_MID, 0.2),
        )
        fig.add_trace(
            go.Scatter(
                x=xi[s], y=nomod[s], line={"color": _PUBU_MID, "width": 2}, name="Unmodified"
            )
        )

        curve_max = max(
            float(np.nanmax(mod + mod_sem)),
            float(np.nanmax(nomod + nomod_sem)),
        )
        y_top = curve_max * 1.1
        pareto_upper = pareto + pareto_sem
        _add_band(fig, xi, pareto_upper, np.full_like(xi, y_top), "rgba(200, 200, 200, 0.3)")
        fig.add_trace(
            go.Scatter(
                x=xi,
                y=pareto_upper,
                line={"color": "black", "width": 1, "dash": "dash"},
                name="Pareto",
            )
        )
        fig.update_yaxes(range=[0, y_top])
    else:
        _add_band(fig, xi, mod - mod_sem, mod + mod_sem, _rgba(_RDPU_MID, 0.2))
        fig.add_trace(
            go.Scatter(x=xi, y=mod, line={"color": _RDPU_MID, "width": 2}, name="Modified")
        )

    fig.add_vline(x=1.0, line={"color": "grey", "dash": "dash", "width": 1}, opacity=0.7)
    title = _title("SNR vs Read Depth", name)
    fig.update_layout(**_base_layout(title, "Relative Total Read Depth", "Mean SNR"))
    fig.update_xaxes(type="log", range=[np.log10(xi[0]), np.log10(xi[-1])])
    return fig
