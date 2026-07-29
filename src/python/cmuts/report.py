"""Generate a self-contained HTML report from a cmuts reactivity HDF5 file.

The report is built purely from the datasets that ``cmuts normalize`` writes; the
HDF5 dataset names are the entire contract between normalize and plot. Each plot
reads the arrays it needs from the file and hands them to a pure plotly builder.
Nothing here imports the normalize pipeline (``ProbingData``, ``Opts``, ...).

Layout (all sections collapsible via ``<details>``):
  * One section per experiment (HDF5 group): summary stats + aggregate plots
    (mutation heatmap, coverage, terminations, SNR-vs-depth); multi-reference
    experiments also get a clipped reactivity heatmap and read-distribution plots.
  * ``embed_all`` (the ``--all`` flag) adds, per experiment, a sequence dropdown
    showing that sequence's profile / coverage / terminations (and MI / correlation
    when present), plus two top-level sections: an extensible multi-plot overlay
    and a difference plot. The coverage/terminations/matrix figures are pre-built
    and swapped client-side; every reactivity profile (per-sequence, multi-plot,
    and difference) is drawn client-side by one shared JS function from the
    embedded reactivity/error arrays, so the three share exactly one code path.

Static images come from plotly's built-in modebar "Download plot as PNG" button,
so no server-side render engine (kaleido/Chrome) is needed.
"""

import argparse
import json
import os
import re
from datetime import date
from typing import Optional

import h5py
import numpy as np
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs

from cmuts.visualize import plotly as P
from cmuts.visualize._transforms import termination_density

# ===========================================================================
# HDF5 dataset-name contract (must match what cmuts normalize writes)
# ===========================================================================
_REACTIVITY = "reactivity"
_READS = "reads"
_ERROR = "error"
_SNR = "SNR"
_HEATMAP = "heatmap"
_COVERAGE = "coverage"
_TERMINATIONS = "terminations"
_MI = "mutual-information"
_COVARIANCE = "covariance"
_META = "meta"  # top-level group holding shared, experiment-independent data
_SEQUENCE = "sequence"  # meta/sequence: per-reference base tokens
_NAME = "name"  # meta/name: per-reference names (FASTA order)
_SNR_XI = "snr-xi"
_SNR_MOD = "snr-mod"
_SNR_MOD_SEM = "snr-mod_sem"
_SNR_NOMOD = "snr-nomod"
_SNR_NOMOD_SEM = "snr-nomod_sem"
_SNR_PARETO = "snr-pareto"
_SNR_PARETO_SEM = "snr-pareto_sem"

# ===========================================================================
# Figure "kinds": on-screen aspect (CSS class) + pinned PNG export size
# ===========================================================================
# Every figure declares a kind. The kind sets a fixed on-screen aspect ratio (via
# the matching .figure CSS class) and a fixed PNG export size (below), so neither
# the layout nor a downloaded image depends on the browser window. "wide" plots
# also span both columns of the plot grid.
#   wide   -- sequence-position plots (5:2), full width
#   std    -- compact line/heatmap plots (3:2), half width
#   square -- L x L matrix plots (data square; card is 5:4 to fit the colorbar)
_FIG_CONFIG = {"displaylogo": False, "responsive": True}
_EXPORT = {"wide": (1000, 400), "std": (600, 400), "square": (650, 520)}

# ===========================================================================
# Page assets: stylesheet and client-side interactivity
# ===========================================================================
_CSS = """
/* Suppress plotly's "Taking snapshot..." notifier on PNG export. */
.plotly-notifier { display: none !important; }
:root {
  color-scheme: light dark;
  --bg: #ffffff; --surface: #f6f6f9; --ink: #1a1a1e; --muted: #6b6b76;
  --border: #e4e4ea; --accent: #c51b8a;
}
@media (prefers-color-scheme: dark) {
  :root { --bg: #14141a; --surface: #1e1e26; --ink: #ececf0; --muted: #9a9aa6;
          --border: #2c2c36; --accent: #e8479f; }
}
:root[data-theme="light"] { --bg: #ffffff; --surface: #f6f6f9; --ink: #1a1a1e;
  --muted: #6b6b76; --border: #e4e4ea; --accent: #c51b8a; }
:root[data-theme="dark"] { --bg: #14141a; --surface: #1e1e26; --ink: #ececf0;
  --muted: #9a9aa6; --border: #2c2c36; --accent: #e8479f; }
* { box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
       max-width: 1200px; margin: 0 auto; padding: 2rem 1.25rem 4rem; line-height: 1.5;
       color: var(--ink); background: var(--bg); }
header.report { border-bottom: 2px solid var(--accent); padding-bottom: 0.7rem; margin-bottom: 1.6rem; }
header.report h1 { font-size: 1.7rem; margin: 0; letter-spacing: -0.01em; }
header.report .subtitle { color: var(--muted); font-size: 0.9rem; margin: 0.3rem 0 0; }
details { border: 1px solid var(--border); border-radius: 10px; margin: 1rem 0;
          padding: 0.3rem 1.1rem 1rem; background: var(--bg); }
details > summary { cursor: pointer; font-size: 1.15rem; font-weight: 650; padding: 0.6rem 0;
                    list-style: none; }
details > summary::-webkit-details-marker { display: none; }
details > summary::before { content: "\\25B8"; color: var(--accent); display: inline-block;
                            width: 1.1em; transition: transform 0.15s ease; }
details[open] > summary::before { transform: rotate(90deg); }
details.plot-details { border: none; border-radius: 0; background: transparent;
                       margin: 0; padding: 0; }
details.plot-details > summary { font-size: 1.05rem; font-weight: 600; color: var(--ink);
                                 padding: 0.5rem 0; }
/* Plots laid out two-up; wide (sequence-position) plots span the full width. */
.plots { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr));
         gap: 0.1rem 1.2rem; align-items: start; grid-auto-flow: row dense; }
@media (max-width: 980px) { .plots { grid-template-columns: 1fr; } }
details.plot-details.wide { grid-column: 1 / -1; }
h3 { font-size: 1.05rem; font-weight: 650; color: var(--muted); text-transform: uppercase;
     letter-spacing: 0.05em; margin: 1.5rem 0 0.4rem; }
.stats-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(150px, 1fr));
              gap: 0.6rem; margin: 0.9rem 0 1.3rem; }
.stat { background: var(--surface); border: 1px solid var(--border); border-radius: 9px;
        padding: 0.7rem 0.9rem; }
.stat-value { font-size: 1.45rem; font-weight: 600; line-height: 1.15; }
.stat-label { font-size: 0.78rem; color: var(--muted); margin-top: 0.25rem; }
/* Figures keep a white background (so downloaded PNGs stay white); the card frames
   that white panel. A fixed aspect-ratio makes the on-screen shape independent of
   the window; the PNG size is pinned separately via toImageButtonOptions. */
.figure { background: #ffffff; border: 1px solid #e6e6ec; border-radius: 10px;
          padding: 6px 10px; margin: 0.2rem 0 0.8rem; box-shadow: 0 1px 4px rgba(0,0,0,0.12); }
.figure.wide { aspect-ratio: 5 / 2; }
.figure.std { aspect-ratio: 3 / 2; }
/* Wider than tall so the plot's width (not height) shrinks to square the data,
   leaving the colorbar at the plot's full, matched height. */
.figure.square { aspect-ratio: 5 / 4; max-width: 760px; margin-inline: auto; }
.figure > *, .figure .plot { width: 100%; height: 100%; }
.plot { width: 100%; }
select, button { margin: 2px 6px 2px 0; padding: 5px 10px; font-size: 0.9rem; color: var(--ink);
                 background: var(--surface); border: 1px solid var(--border); border-radius: 6px; }
button { cursor: pointer; }
.row { margin: 5px 0; }
.seqview-wrap { margin: 4px 0 14px; }
.seqview-label { font-size: 0.78rem; font-weight: 600; color: var(--muted);
                 text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 4px; }
.seqview { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12px;
           line-height: 1.6; word-break: break-all; border: 1px solid #2c2c36;
           border-radius: 6px; padding: 6px 8px; background: #14141a; }
.seqview .base { display: inline-block; width: 1.15ch; text-align: center; font-weight: 600; }
#diff-rows { display: flex; align-items: center; flex-wrap: wrap; gap: 8px; }
#diff-rows .row { margin: 0; }
#diff-rows select { margin: 0; }
.diff-op { color: var(--muted); }
.structure pre { white-space: pre-wrap; background: var(--surface); border: 1px solid var(--border);
                 padding: 0.75rem; border-radius: 8px; overflow-x: auto; }
""".strip()

# Client-side interactivity. Reads three embedded JSON blobs:
#   CMUTS_FIGS  -- pre-built per-sequence figures, {slug: {ref: {key: figure}}}
#   CMUTS_DATA  -- raw arrays,   {slug: {label, reactivity: [[..]], error: [[..]], n}}
#   CMUTS_EXPS  -- [{slug, label, n}] in report order
# The per-sequence dropdown swaps pre-built figures; the multi-plot and
# difference plot compose reactivity line traces from CMUTS_DATA. The EXPORT and
# PSEQ_KIND maps below mirror the Python _EXPORT / _PER_SEQ tables.
_REPORT_JS = r"""
(function () {
  var COLORS = ["#c51b8a", "#2c7fb8", "#31a354", "#e6550d", "#756bb1", "#636363"];
  // PNG export size per figure kind, so downloads don't depend on the window.
  var EXPORT = {wide: [1000, 400], std: [600, 400], square: [650, 520]};
  var PSEQ_KIND = {profile: "wide", coverage: "std", terminations: "std",
                   mi: "square", correlation: "square"};
  function slugify(s) {
    return String(s).replace(/[^A-Za-z0-9_-]+/g, "-").replace(/^-+|-+$/g, "") || "plot";
  }
  function cfgFor(kind, title) {
    var e = EXPORT[kind] || EXPORT.std;
    var opts = {format: "png", scale: 2, width: e[0], height: e[1]};
    if (title) opts.filename = slugify(title);
    return {displaylogo: false, responsive: true, toImageButtonOptions: opts};
  }

  function readJSON(id) {
    var el = document.getElementById(id);
    return el ? JSON.parse(el.textContent) : null;
  }
  var FIGS = readJSON("cmuts-figs") || {};
  var DATA = readJSON("cmuts-data") || {};
  var EXPS = readJSON("cmuts-exps") || [];
  var SEQS = readJSON("cmuts-seq") || [];   // per-reference base strings, FASTA order
  var SINGLE_EXP = EXPS.length === 1;

  // --- colored sequence view: each base's glyph shaded white -> purple ---------
  function seqColor(v) {
    if (v === null || v === undefined || isNaN(v)) return "var(--muted)";
    var t = Math.max(0, Math.min(1, v));   // #762a83 at t=1, white at t=0
    return "rgb(" + Math.round(255 + (118 - 255) * t) + ","
                  + Math.round(255 + (42 - 255) * t) + ","
                  + Math.round(255 + (131 - 255) * t) + ")";
  }
  function renderSeqView(div, letters, react) {
    if (!div) return;
    if (!letters) { div.innerHTML = ""; return; }
    var html = "";
    for (var i = 0; i < letters.length; i++) {
      html += '<span class="base" style="color:' + seqColor(react ? react[i] : null) + '">'
            + letters[i] + "</span>";
    }
    div.innerHTML = html;
  }
  // Single-reference sequence views have no dropdown; fill them once on load.
  function bindStaticSeqviews() {
    var els = document.querySelectorAll(".js-seqview");
    for (var i = 0; i < els.length; i++) {
      var slug = els[i].getAttribute("data-group");
      var ref = parseInt(els[i].getAttribute("data-ref") || "0", 10);
      var d = DATA[slug];
      renderSeqView(els[i], SEQS[ref], d ? d.reactivity[ref] : null);
    }
  }

  function baseLayout(title, ylabel) {
    return {title: {text: title}, xaxis: {title: {text: "Residue"}},
            yaxis: {title: {text: ylabel}}, margin: {t: 40, r: 20, b: 45, l: 55}};
  }
  function xrange(n) { var x = new Array(n); for (var i = 0; i < n; i++) x[i] = i; return x; }

  // --- shared profile renderer: line + surround error band ------------------
  // Used for every reactivity profile: per-sequence, multi-plot overlay, and the
  // difference plot. The band is a single self-closing polygon (upper then lower
  // reversed) so it composes correctly when several profiles are overlaid.
  var PROFILE_COLOR = "#c51b8a";
  function hexToRgba(hex, a) {
    var n = parseInt(hex.slice(1), 16);
    return "rgba(" + ((n >> 16) & 255) + "," + ((n >> 8) & 255) + "," + (n & 255) + "," + a + ")";
  }
  function defined(v) { return v !== null && v !== undefined && !isNaN(v); }
  function profileTraces(y, err, name, color, errorLabel) {
    color = color || PROFILE_COLOR;
    var x = xrange(y.length);
    var traces = [];
    if (err) {
      // One filled polygon per contiguous run of defined points, so the band
      // breaks at gaps (masked positions) instead of filling across them. Only
      // the first run carries the legend entry; labelling is used for a single
      // profile (errorLabel given) -- overlaid bands share their line's colour,
      // so labelling each would only clutter the legend.
      var shown = false, i = 0;
      while (i < y.length) {
        if (!defined(y[i])) { i++; continue; }
        var s = i;
        while (i < y.length && defined(y[i])) i++;
        var bx = [], by = [];
        for (var k = s; k < i; k++) { bx.push(x[k]); by.push(y[k] + (err[k] || 0)); }
        for (var j = i - 1; j >= s; j--) { bx.push(x[j]); by.push(y[j] - (err[j] || 0)); }
        traces.push({x: bx, y: by, fill: "toself", fillcolor: hexToRgba(color, 0.2),
                     line: {width: 0}, mode: "lines",
                     name: errorLabel || "", showlegend: !!errorLabel && !shown,
                     hoverinfo: "skip"});
        shown = true;
      }
    }
    traces.push({x: x, y: y, mode: "lines", line: {color: color, width: 1.5},
                 name: name || "Reactivity"});
    return traces;
  }
  // A lone profile: legend labels the line and its standard-error band; the title
  // (which names the experiment/reference) also seeds the download filename.
  function renderProfile(div, y, err, title, color) {
    var layout = baseLayout(title, "Reactivity");
    layout.showlegend = true;
    // Inside the top-right corner (the profile is near-zero there), so the legend
    // overlays the plot instead of reserving a white column beside it.
    // Reversed so the line ("Reactivity") lists above its band ("Std. error"),
    // while the band is still drawn first (behind the line).
    layout.legend = {x: 0.99, xanchor: "right", y: 0.99, yanchor: "top",
                     bgcolor: "rgba(0,0,0,0)", borderwidth: 0, traceorder: "reversed"};
    Plotly.react(div, profileTraces(y, err, "Reactivity", color, "Std. error"),
                 layout, cfgFor("wide", title));
  }

  // --- per-experiment sequence dropdown: swap pre-built figures --------------
  function bindDropdowns() {
    var sels = document.querySelectorAll(".seqselect");
    for (var k = 0; k < sels.length; k++) {
      (function (sel) {
        var slug = sel.getAttribute("data-group");
        function render() {
          var d = DATA[slug]; if (!d) return;
          // Per-sequence stat tiles for the selected reference.
          if (d.seqstats && d.seqstats[sel.value]) {
            var ss = d.seqstats[sel.value];
            for (var s in ss) {
              var st = document.getElementById("pseqstat-" + slug + "-" + s);
              if (st) st.textContent = ss[s];
            }
          }
          // Profile: built client-side with the shared renderer.
          var pdiv = document.getElementById("perseq-" + slug + "-profile");
          if (pdiv) renderProfile(pdiv, d.reactivity[sel.value], d.error[sel.value],
                                  "Reactivity profile (" + d.label + " · "
                                  + refName(slug, sel.value) + ")", PROFILE_COLOR);
          // Colored sequence for the selected reference.
          renderSeqView(document.getElementById("perseq-" + slug + "-seq"),
                        SEQS[sel.value], d.reactivity[sel.value]);
          // Remaining per-sequence figures: pre-built figure swaps.
          var byRef = FIGS[slug];
          var fig = byRef ? byRef[sel.value] : null;
          if (fig) {
            for (var key in fig) {
              var div = document.getElementById("perseq-" + slug + "-" + key);
              var fttl = fig[key].layout && fig[key].layout.title ? fig[key].layout.title.text : key;
              if (div) Plotly.react(div, fig[key].data, fig[key].layout,
                                    cfgFor(PSEQ_KIND[key] || "std", fttl));
            }
          }
        }
        sel.addEventListener("change", render);
        render();
      })(sels[k]);
    }
  }

  // --- experiment -> sequence chained selects --------------------------------
  function refName(slug, i) {
    var d = DATA[slug];
    return (d && d.names && d.names[i]) ? d.names[i] : ("Reference " + (i + 1));
  }
  function fillSequences(seqSel, slug) {
    var n = DATA[slug] ? DATA[slug].n : 0;
    var opts = "";
    for (var i = 0; i < n; i++) opts += '<option value="' + i + '">' + refName(slug, i) + "</option>";
    seqSel.innerHTML = opts;
  }
  function makeRow(onChange, removable) {
    var row = document.createElement("div");
    row.className = "row";
    // Experiment selector only when there is more than one experiment.
    var expSel = null;
    if (!SINGLE_EXP) {
      expSel = document.createElement("select");
      for (var i = 0; i < EXPS.length; i++)
        expSel.innerHTML += '<option value="' + EXPS[i].slug + '">' + EXPS[i].label + "</option>";
      row.appendChild(expSel);
    }
    var seqSel = document.createElement("select");
    row.appendChild(seqSel);
    function slug() { return SINGLE_EXP ? EXPS[0].slug : expSel.value; }
    // The sequence selector is only useful when the experiment has >1 reference.
    function refreshSeq() {
      var s = slug(); var n = DATA[s] ? DATA[s].n : 0;
      if (n > 1) {
        fillSequences(seqSel, s); seqSel.style.display = "";
      } else {
        seqSel.innerHTML = '<option value="0"></option>'; seqSel.value = "0";
        seqSel.style.display = "none";
      }
    }
    if (expSel) expSel.addEventListener("change", function () { refreshSeq(); onChange(); });
    seqSel.addEventListener("change", onChange);
    refreshSeq();
    if (removable) {
      var rm = document.createElement("button");
      rm.textContent = "−";
      rm.addEventListener("click", function () { row.parentNode.removeChild(row); onChange(); });
      row.appendChild(rm);
    }
    row._sel = function () { return {slug: slug(), ref: parseInt(seqSel.value || "0", 10)}; };
    return row;
  }
  // Legend label mirrors the dropdowns: the experiment name only when there is
  // more than one experiment, the FASTA name only for multi-reference experiments.
  function traceLabel(sel) {
    var parts = [];
    if (!SINGLE_EXP) parts.push(DATA[sel.slug].label);
    if (DATA[sel.slug].n > 1) parts.push(refName(sel.slug, sel.ref));
    return parts.join(" · ");
  }
  // --- multi-plot: arbitrarily many overlaid reactivity profiles -------------
  function bindMultiPlot() {
    var rowsEl = document.getElementById("multiplot-rows");
    var addBtn = document.getElementById("multiplot-add");
    var plot = document.getElementById("multiplot");
    if (!rowsEl || !plot) return;
    function redraw() {
      var traces = [];
      var rows = rowsEl.querySelectorAll(".row");
      for (var i = 0; i < rows.length; i++) {
        var sel = rows[i]._sel(), d = DATA[sel.slug];
        traces = traces.concat(profileTraces(
          d.reactivity[sel.ref], d.error[sel.ref], traceLabel(sel), COLORS[i % COLORS.length]));
      }
      var layout = baseLayout("Reactivity overlay", "Reactivity");
      layout.showlegend = true;
      layout.legend = {x: 0.99, xanchor: "right", y: 0.99, yanchor: "top",
                       bgcolor: "rgba(0,0,0,0)", borderwidth: 0};
      Plotly.react(plot, traces, layout, cfgFor("wide", "Reactivity overlay"));
    }
    function addRow() { rowsEl.appendChild(makeRow(redraw, true)); redraw(); }
    addBtn.addEventListener("click", addRow);
    addRow();  // start with one row
  }

  // --- difference plot: exactly two sequences, plot A - B --------------------
  function bindDiffPlot() {
    var rowsEl = document.getElementById("diff-rows");
    var plot = document.getElementById("diffplot");
    if (!rowsEl || !plot) return;
    var rowA = makeRow(redraw, false), rowB = makeRow(redraw, false);
    var op = document.createElement("span");
    op.className = "diff-op"; op.textContent = "minus";
    rowsEl.appendChild(rowA); rowsEl.appendChild(op); rowsEl.appendChild(rowB);
    function redraw() {
      var a = rowA._sel(), b = rowB._sel();
      var da = DATA[a.slug], db = DATA[b.slug];
      var ya = da.reactivity[a.ref], yb = db.reactivity[b.ref];
      var ea = da.error[a.ref], eb = db.error[b.ref];
      var n = Math.min(ya.length, yb.length);
      var diff = [], err = [];
      for (var i = 0; i < n; i++) {
        // A masked position in either profile leaves the difference undefined;
        // push NaN so the line/band break there (embedded gaps arrive as null,
        // which would otherwise coerce to 0 in the subtraction).
        if (!defined(ya[i]) || !defined(yb[i])) { diff.push(NaN); err.push(NaN); continue; }
        diff.push(ya[i] - yb[i]);
        err.push(Math.sqrt((ea[i] || 0) * (ea[i] || 0) + (eb[i] || 0) * (eb[i] || 0)));
      }
      var label = da.label + " · " + refName(a.slug, a.ref) + " − "
                + db.label + " · " + refName(b.slug, b.ref);
      // A single line, so no legend (the label would only be long clutter).
      var layout = baseLayout("Difference", "Reactivity Difference");
      layout.showlegend = false;
      Plotly.react(plot, profileTraces(diff, err, label), layout,
                   cfgFor("wide", "Difference " + label));
    }
    redraw();
  }

  // Sections are collapsed by default; a plot drawn while hidden has zero size,
  // so resize any plots inside a <details> the first time it is opened.
  function bindCollapsibleResize() {
    var ds = document.querySelectorAll("details");
    for (var i = 0; i < ds.length; i++) {
      ds[i].addEventListener("toggle", function () {
        if (!this.open) return;
        var plots = this.querySelectorAll(".plotly-graph-div");
        for (var j = 0; j < plots.length; j++) {
          if (plots[j].offsetWidth > 0) { try { Plotly.Plots.resize(plots[j]); } catch (e) {} }
        }
      });
    }
  }

  // Single-reference profiles have no dropdown; render them once on load.
  function bindStaticProfiles() {
    var els = document.querySelectorAll(".js-profile");
    for (var i = 0; i < els.length; i++) {
      var slug = els[i].getAttribute("data-group");
      var ref = parseInt(els[i].getAttribute("data-ref") || "0", 10);
      var d = DATA[slug];
      if (d) renderProfile(els[i], d.reactivity[ref], d.error[ref],
                           "Reactivity profile (" + (d.label || slug) + ")", PROFILE_COLOR);
    }
  }

  document.addEventListener("DOMContentLoaded", function () {
    bindStaticProfiles(); bindStaticSeqviews(); bindDropdowns(); bindMultiPlot(); bindDiffPlot();
    bindCollapsibleResize();
  });
})();
""".strip()


# ===========================================================================
# HDF5 read helpers
# ===========================================================================
def _slug(name: str) -> str:
    """HTML-id-safe token for a group name (used in element ids / data-group)."""
    return re.sub(r"[^A-Za-z0-9_-]", "-", name) or "root"


def _arr(grp: h5py.Group, name: str) -> Optional[np.ndarray]:
    """The named dataset as an array, or None if the group lacks it."""
    return np.asarray(grp[name]) if name in grp else None


def _str_list(container, name: str) -> Optional[list]:
    """A dataset of (byte)strings decoded to ``list[str]``, or None if absent."""
    if name not in container:
        return None
    return [s.decode() if isinstance(s, bytes) else str(s) for s in container[name]]


def _meta_arr(f: h5py.File, name: str) -> Optional[np.ndarray]:
    """A shared dataset from the ``meta`` group, falling back to a legacy top-level
    dataset of the same name (files written before the meta group existed)."""
    if _META in f and name in f[_META]:
        return np.asarray(f[_META][name])
    if name in f and not isinstance(f[name], h5py.Group):
        return np.asarray(f[name])
    return None


def _names_list(f: h5py.File) -> Optional[list]:
    """Per-reference names from ``meta/name`` (or the legacy top-level ``names``)."""
    if _META in f and _NAME in f[_META]:
        return _str_list(f[_META], _NAME)
    return _str_list(f, "names")


# Inverse of the default base->token map (see internal._BASE_TOKEN): 0/1/2/3 -> ACGU.
_TOKEN_BASE = {0: "A", 1: "C", 2: "G", 3: "U"}


def _decode_sequences(f: h5py.File) -> Optional[list]:
    """The ``meta/sequence`` token array decoded to per-reference base strings (FASTA
    order), or None if absent. Trailing padding (-1) is dropped; any other non-base
    token renders as ``N``. Powers the colored per-base sequence view."""
    tokens = _meta_arr(f, _SEQUENCE)
    if tokens is None:
        return None
    seqs = []
    for row in np.atleast_2d(tokens):
        end = len(row)
        while end > 0 and int(row[end - 1]) < 0:  # trim trailing padding
            end -= 1
        seqs.append("".join(_TOKEN_BASE.get(int(t), "N") for t in row[:end]))
    return seqs


def _ref_label(names: Optional[list], i: int) -> str:
    """Label reference ``i`` by its FASTA name, falling back to a 1-based index."""
    if names and i < len(names) and names[i]:
        return str(names[i])
    return f"Reference {i + 1}"


def _first_ref(arr: np.ndarray) -> np.ndarray:
    """First reference's row of a per-reference array; tolerates legacy 1-D."""
    return arr[0] if arr.ndim == 2 else arr


def _coverage_aggregate(coverage: np.ndarray) -> np.ndarray:
    """Mean coverage over references. Tolerates legacy 1-D (already aggregate)."""
    return coverage.mean(0) if coverage.ndim == 2 else coverage


def _clean_2d(arr: np.ndarray) -> list:
    """2-D float list with NaN -> None, so it survives JSON.parse in the browser."""
    a = np.asarray(arr, dtype=float)
    return [[None if np.isnan(v) else float(v) for v in row] for row in a]


# ===========================================================================
# Summary statistics
# ===========================================================================
def _stats(grp: h5py.Group) -> list:
    """Experiment-wide summary statistics as (label, formatted value) rows."""
    reactivity = np.asarray(grp[_REACTIVITY])
    reads = np.asarray(grp[_READS])
    error = _arr(grp, _ERROR)
    snr = _arr(grp, _SNR)

    n_refs = reactivity.shape[0]
    seq_len = reactivity.shape[1] if reactivity.ndim > 1 else 0
    valid = np.isfinite(reactivity)

    rows = [
        ("Reference" if n_refs == 1 else "References", f"{n_refs:,}"),
        ("Reference length", f"{seq_len:,}"),
        ("Total reads", f"{int(reads.sum()):,}"),
    ]
    if n_refs > 1:
        # Per-reference mean/median are redundant with total for a single reference.
        rows.append(("Mean reads per reference", f"{np.mean(reads):,.1f}"))
        rows.append(("Median reads per reference", f"{int(np.median(reads)):,}"))
    if valid.any():
        rows.append(("Mean reactivity", f"{np.mean(reactivity[valid]):.3f}"))
        if error is not None:
            rows.append(("Mean error", f"{np.mean(error[valid]):.3f}"))
        if snr is not None:
            # SNR is per reference: its mean over one reference is just that SNR,
            # and the "fraction with SNR > 1" is only meaningful across references.
            rows.append(("SNR" if n_refs == 1 else "Mean SNR", f"{np.mean(snr):.2f}"))
            if n_refs > 1:
                rows.append(("SNR > 1", f"{np.mean(snr > 1):.1%}"))
    dropout = float(np.mean(reads == 0))
    if dropout > 0:
        rows.append(("Dropout fraction", f"{dropout:.1%}"))
    return rows


# Per-sequence stat tiles: (embed key, label). Values are computed by _seq_stats
# and written into the tiles client-side when a reference is selected.
_PER_SEQ_STATS = [
    ("reads", "Reads"),
    ("snr", "SNR"),
    ("reactivity", "Mean reactivity"),
    ("error", "Mean error"),
]


def _seq_stats(grp: h5py.Group, n_refs: int) -> list:
    """Per-reference stat values (formatted strings) keyed by ``_PER_SEQ_STATS``."""
    reads = np.asarray(grp[_READS])
    snr = _arr(grp, _SNR)
    reactivity = np.asarray(grp[_REACTIVITY])
    error = _arr(grp, _ERROR)
    out = []
    for i in range(n_refs):
        row = reactivity[i]
        valid = np.isfinite(row)
        out.append(
            {
                "reads": f"{int(reads[i]):,}",
                "snr": f"{float(snr[i]):.2f}" if snr is not None else "—",
                "reactivity": f"{float(np.nanmean(row)):.3f}" if valid.any() else "—",
                "error": f"{float(np.nanmean(error[i])):.3f}"
                if (error is not None and valid.any())
                else "—",
            }
        )
    return out


# ===========================================================================
# HTML fragment helpers
# ===========================================================================
def _stats_grid_html(rows: list) -> str:
    """Summary statistics as a responsive grid of stat tiles (value over label)."""
    tiles = "".join(
        f'<div class="stat"><div class="stat-value">{value}</div>'
        f'<div class="stat-label">{label}</div></div>'
        for label, value in rows
    )
    return f'<div class="stats-grid">{tiles}</div>'


def _per_seq_stats_grid_html(slug: str) -> str:
    """Empty per-sequence stat tiles; values filled in by the dropdown JS."""
    tiles = "".join(
        f'<div class="stat"><div class="stat-value" id="pseqstat-{slug}-{key}"></div>'
        f'<div class="stat-label">{label}</div></div>'
        for key, label in _PER_SEQ_STATS
    )
    return f'<div class="stats-grid">{tiles}</div>'


def _fig_html(fig: go.Figure, kind: str) -> str:
    """Inline a rendered figure with a pinned PNG export size for its kind, and a
    descriptive download filename taken from the figure's title."""
    w, h = _EXPORT.get(kind, _EXPORT["std"])
    filename = re.sub(r"-+", "-", _slug(fig.layout.title.text or "plot")).strip("-") or "plot"
    cfg = {
        **_FIG_CONFIG,
        "toImageButtonOptions": {
            "format": "png",
            "scale": 2,
            "width": w,
            "height": h,
            "filename": filename,
        },
    }
    html = fig.to_html(full_html=False, include_plotlyjs=False, config=cfg)
    return f'<div class="figure {kind}">{html}</div>'


def _plot_block(title: str, inner_html: str, kind: str) -> str:
    """A single plot as its own collapsible block, toggled by its (large) title.

    A ``wide`` kind makes the block span both columns of the plot grid.
    """
    cls = "plot-details wide" if kind == "wide" else "plot-details"
    return f'<details class="{cls}"><summary>{title}</summary>{inner_html}</details>'


def _figure_block(title: str, fig: go.Figure, kind: str) -> str:
    """Collapsible block wrapping a fully rendered figure."""
    return _plot_block(title, _fig_html(fig, kind), kind)


def _placeholder_block(title: str, div_id: str, kind: str) -> str:
    """Collapsible block with an empty plot div, populated client-side."""
    inner = f'<div class="figure {kind}"><div class="plot" id="{div_id}"></div></div>'
    return _plot_block(title, inner, kind)


def _js_profile_block(title: str, slug: str, ref: int) -> str:
    """Reactivity-profile block rendered client-side by the shared profile
    function (from embedded arrays), so it matches the multi-plot / difference
    plot. ``.js-profile`` divs are drawn once on load; the per-sequence dropdown
    targets ``perseq-<slug>-profile`` instead (see the report JS)."""
    inner = (
        f'<div class="figure wide">'
        f'<div class="plot js-profile" data-group="{slug}" data-ref="{ref}"></div></div>'
    )
    return _plot_block(title, inner, "wide")


def _seqview_html(slug: str, ref: Optional[int]) -> str:
    """The per-base colored sequence strip (white -> purple by reactivity).

    ``ref`` set (single reference) -> a ``.js-seqview`` filled once on load; ``ref``
    None (multi-reference) -> an id'd div the dropdown refills for the selection."""
    inner = (
        f'<div class="seqview js-seqview" data-group="{slug}" data-ref="{ref}"></div>'
        if ref is not None
        else f'<div class="seqview" id="perseq-{slug}-seq"></div>'
    )
    return f'<div class="seqview-wrap"><div class="seqview-label">Sequence</div>{inner}</div>'


def _plots_grid(blocks: list) -> str:
    """Wrap plot blocks in the two-column plot grid."""
    return '<div class="plots">' + "".join(blocks) + "</div>"


def _section(title: str, body: str) -> str:
    """A top-level collapsible section (experiment / comparison / structure)."""
    return f"<details><summary>{title}</summary>{body}</details>"


# ===========================================================================
# Figure builders (HDF5 group -> plotly figures, tagged with a kind)
# ===========================================================================
def _snr_scaling_fig(grp: h5py.Group, name: str) -> Optional[go.Figure]:
    """SNR-vs-read-depth figure, or None if the curves are absent."""
    xi = _arr(grp, _SNR_XI)
    if xi is None or _SNR_MOD not in grp or _SNR_MOD_SEM not in grp:
        return None
    return P.plot_snr_scaling(
        xi,
        np.asarray(grp[_SNR_MOD]),
        np.asarray(grp[_SNR_MOD_SEM]),
        _arr(grp, _SNR_NOMOD),
        _arr(grp, _SNR_NOMOD_SEM),
        _arr(grp, _SNR_PARETO),
        _arr(grp, _SNR_PARETO_SEM),
        name,
    )


def _aggregate_figures(grp: h5py.Group, name: str) -> list:
    """Reference-independent figures as (title, figure, kind); skips missing data."""
    reads = _arr(grp, _READS)
    n_refs = int(reads.shape[0]) if reads is not None else 0
    figs = []

    heatmap = _arr(grp, _HEATMAP)
    if heatmap is not None:
        figs.append(("Mutation heatmap", P.plot_heatmap(heatmap, name), "std"))

    coverage = _arr(grp, _COVERAGE)
    if coverage is not None and reads is not None:
        figs.append(
            ("Coverage", P.plot_coverage(_coverage_aggregate(coverage), reads, name), "std")
        )

    terminations = _arr(grp, _TERMINATIONS)
    if terminations is not None:
        figs.append(
            ("Terminations", P.plot_termination(termination_density(terminations), name), "std")
        )

    snr = _snr_scaling_fig(grp, name)
    if snr is not None:
        figs.append(("SNR vs read depth", snr, "std"))

    if n_refs > 1:  # read distributions and the profile heatmap need many references
        assert reads is not None  # n_refs > 1 implies reads was present
        reactivity = _arr(grp, _REACTIVITY)
        if reactivity is not None:
            figs.append(("Reactivity heatmap", P.plot_reactivity_heatmap(reactivity, name), "wide"))
        figs.append(("Reads per reference", P.plot_read_hist(reads, name), "std"))
        figs.append(("Reads per block", P.plot_reads_per_block(reads, name), "std"))
    return figs


def _single_reference_figures(grp: h5py.Group, name: str) -> list:
    """Server-rendered (title, figure, kind) plots for a one-reference experiment.

    The reactivity profile is rendered client-side (shared with the multi-plot /
    difference plot), so it is not included here. With one reference the aggregate
    and per-sequence views coincide, so every plot is shown once (and without
    ``embed_all``, since it is bounded).
    """
    reads = np.asarray(grp[_READS])
    figs = []
    heatmap = _arr(grp, _HEATMAP)
    if heatmap is not None:
        figs.append(("Mutation heatmap", P.plot_heatmap(heatmap, name), "std"))
    coverage = _arr(grp, _COVERAGE)
    if coverage is not None:
        figs.append(("Coverage", P.plot_coverage(_first_ref(coverage), reads, name), "std"))
    terminations = _arr(grp, _TERMINATIONS)
    if terminations is not None:
        figs.append(("Terminations", P.plot_termination(_first_ref(terminations), name), "std"))
    snr = _snr_scaling_fig(grp, name)
    if snr is not None:
        figs.append(("SNR vs read depth", snr, "std"))
    mi = _arr(grp, _MI)
    if mi is not None:
        figs.append(("Mutual information", P.plot_mi(mi[0], name), "square"))
    covariance = _arr(grp, _COVARIANCE)
    if covariance is not None:
        figs.append(("Correlation", P.plot_correlation(covariance[0], name), "square"))
    return figs


# Per-sequence plots: embed key -> (section title, figure aspect kind).
# The aspect kinds mirror the JS PSEQ_KIND map in _REPORT_JS.
_PER_SEQ = {
    "profile": ("Reactivity profile", "wide"),
    "coverage": ("Coverage", "std"),
    "terminations": ("Terminations", "std"),
    "mi": ("Mutual information", "square"),
    "correlation": ("Correlation", "square"),
}


def _per_sequence_figures(grp: h5py.Group, name: str, ref: int) -> dict:
    """Pre-built figures for one reference, keyed by ``_PER_SEQ``; only ones with
    data. The profile is excluded -- it is rendered client-side by the shared
    profile function (see the report JS)."""
    reads = np.asarray(grp[_READS])
    figs = {}
    coverage = _arr(grp, _COVERAGE)
    if coverage is not None and coverage.ndim == 2:
        figs["coverage"] = P.plot_coverage(coverage[ref], reads[ref], name)
    terminations = _arr(grp, _TERMINATIONS)
    if terminations is not None:
        figs["terminations"] = P.plot_termination(terminations[ref], name)
    mi = _arr(grp, _MI)
    if mi is not None:
        figs["mi"] = P.plot_mi(mi[ref], name)
    covariance = _arr(grp, _COVARIANCE)
    if covariance is not None:
        figs["correlation"] = P.plot_correlation(covariance[ref], name)
    return figs


# ===========================================================================
# Per-experiment sections
# ===========================================================================
def _single_reference_section(grp: h5py.Group, name: str, slug: str, has_seq: bool) -> str:
    """Stats + colored sequence + a flat grid of plots for a one-reference
    experiment. The profile is client-rendered (from the embedded arrays); the rest
    are server-rendered."""
    blocks = [_js_profile_block("Reactivity profile", slug, 0)]
    blocks += [_figure_block(t, fig, kind) for t, fig, kind in _single_reference_figures(grp, name)]
    seq = _seqview_html(slug, 0) if has_seq else ""
    return _stats_grid_html(_stats(grp)) + seq + _plots_grid(blocks)


def _aggregate_section(grp: h5py.Group, name: str) -> str:
    """ "Aggregate data" heading + stats + aggregate plots (multi-reference)."""
    blocks = [_figure_block(t, fig, kind) for t, fig, kind in _aggregate_figures(grp, name)]
    return "<h3>Aggregate data</h3>" + _stats_grid_html(_stats(grp)) + _plots_grid(blocks)


def _per_sequence_section(
    grp: h5py.Group,
    name: str,
    slug: str,
    names: Optional[list],
    n_refs: int,
    has_seq: bool,
) -> tuple:
    """ "Per-sequence data" section (multi-reference, ``embed_all``).

    Returns (html, figs_by_ref) where figs_by_ref maps each reference index to its
    pre-built figures (as plotly JSON) for the dropdown to swap in client-side. The
    profile block is client-rendered; the matrix blocks are pre-built swaps.
    """
    by_ref = {}
    keys = []
    for i in range(n_refs):
        figs = _per_sequence_figures(grp, name, i)
        by_ref[i] = {k: json.loads(v.to_json()) for k, v in figs.items()}
        keys = list(figs.keys())

    options = "".join(f'<option value="{i}">{_ref_label(names, i)}</option>' for i in range(n_refs))
    blocks = [
        _placeholder_block(_PER_SEQ["profile"][0], f"perseq-{slug}-profile", _PER_SEQ["profile"][1])
    ]
    blocks += [
        _placeholder_block(_PER_SEQ[k][0], f"perseq-{slug}-{k}", _PER_SEQ[k][1]) for k in keys
    ]
    html = (
        f'<h3>Per-sequence data</h3><select class="seqselect" data-group="{slug}">{options}</select>'
        + _per_seq_stats_grid_html(slug)
        + (_seqview_html(slug, None) if has_seq else "")
        + _plots_grid(blocks)
    )
    return html, by_ref


def _embed_data(grp: h5py.Group, label: str, names: Optional[list], n_refs: int) -> dict:
    """Raw per-reference arrays + labels + stats for the multi-plot / diff JS."""
    reactivity = np.asarray(grp[_REACTIVITY])
    error = _arr(grp, _ERROR)
    return {
        "label": label,
        "reactivity": _clean_2d(reactivity),
        "error": _clean_2d(error) if error is not None else _clean_2d(reactivity * 0),
        "n": int(n_refs),
        "names": [_ref_label(names, i) for i in range(n_refs)],
        "seqstats": _seq_stats(grp, n_refs),
    }


# ===========================================================================
# Cross-experiment sections and document assembly
# ===========================================================================
def _can_compare(embed_all: bool, exps: list) -> bool:
    """Whether the multi-plot / difference sections have >=2 comparable profiles:
    more than one experiment, or a single experiment with multiple references."""
    return bool(embed_all and exps and (len(exps) > 1 or exps[0]["n"] > 1))


def _comparison_sections() -> list:
    """The multi-plot and difference sections (static shells; JS wires them up)."""
    return [
        _section(
            "Multi-plot",
            '<div id="multiplot-rows"></div>'
            '<button id="multiplot-add">+ Add sequence</button>'
            '<div class="figure wide"><div class="plot" id="multiplot"></div></div>',
        ),
        _section(
            "Difference plot",
            '<div id="diff-rows"></div>'
            '<div class="figure wide"><div class="plot" id="diffplot"></div></div>',
        ),
    ]


def _embed_tags(figs_embed: dict, data_embed: dict, exps: list, seq_embed: Optional[list]) -> str:
    """The ``application/json`` blobs the client JS reads on load (only the non-empty
    ones): pre-built per-sequence figures, raw arrays, experiment list, and the
    shared per-reference base strings (indexed in FASTA order)."""

    def tag(name: str, obj: object) -> str:
        return f'<script type="application/json" id="{name}">{json.dumps(obj)}</script>'

    parts = []
    if figs_embed:
        parts.append(tag("cmuts-figs", figs_embed))
    if data_embed:
        parts.append(tag("cmuts-data", data_embed))
    if exps:
        parts.append(tag("cmuts-exps", exps))
    if seq_embed:
        parts.append(tag("cmuts-seq", seq_embed))
    return ("\n".join(parts) + "\n") if parts else ""


def _document(sections_html: str, embed_tags: str, n_experiments: int) -> str:
    """Assemble the complete HTML document around the section bodies."""
    today = date.today()
    date_str = f"{today:%B} {today.day}, {today.year}"
    plural = "s" if n_experiments != 1 else ""
    subtitle = f"Generated {date_str} · {n_experiments} experiment{plural}"
    return (
        "<!doctype html>\n"
        '<html><head><meta charset="utf-8">'
        f"<title>cmuts report {date_str}</title>\n"
        f"<style>{_CSS}</style>\n"
        f"<script>{get_plotlyjs()}</script>\n"
        "</head><body>\n"
        '<header class="report"><h1>cmuts report</h1>'
        f'<p class="subtitle">{subtitle}</p></header>\n'
        + sections_html
        + "\n"
        + embed_tags
        + f"<script>{_REPORT_JS}</script>\n"
        "</body></html>\n"
    )


# ===========================================================================
# Public API
# ===========================================================================
def build(
    h5_path: str,
    *,
    embed_all: bool = False,
    group: Optional[str] = None,
    structure: Optional[str] = None,
) -> str:
    """Build a self-contained HTML report from a reactivity HDF5 file.

    Args:
        h5_path: Path to a reactivity HDF5 written by ``cmuts normalize``.
        embed_all: Also embed every reference's per-sequence figures (dropdown) and
            the multi-plot / difference sections. Off by default so the report
            stays small on reference-heavy runs.
        group: Restrict to a single group (default: every top-level group).
        structure: Optional pre-rendered markdown/text block (e.g. ChimeraX
            commands) appended as a "Structure visualization" section.

    Returns:
        A complete HTML document as a string.
    """
    sections = []
    figs_embed = {}  # slug -> {ref -> {key -> plotly figure json}}   (per-sequence dropdown)
    data_embed = {}  # slug -> raw arrays + labels + stats            (multi-plot / diff)
    exps = []  # [{slug, label, n}] in report order

    with h5py.File(h5_path, "r") as f:
        if group is not None:
            groups = [group]
        else:
            groups = [k for k in f if isinstance(f[k], h5py.Group) and k != _META] or [""]
        names = _names_list(f)
        sequences = _decode_sequences(f)  # per-reference base strings, shared across groups
        has_seq = sequences is not None

        for g in groups:
            grp = f if g == "" else f[g]
            slug = _slug(g)
            label = g or "Profile"
            n_refs = int(np.asarray(grp[_REACTIVITY]).shape[0])
            single = n_refs == 1

            if single:
                body = _single_reference_section(grp, g, slug, has_seq)
            else:
                body = _aggregate_section(grp, g)
                if embed_all:
                    per_seq_html, figs_embed[slug] = _per_sequence_section(
                        grp, g, slug, names, n_refs, has_seq
                    )
                    body += per_seq_html
            sections.append(_section(label, body))

            # Every rendered profile is client-side, so embed its arrays: always
            # for a single reference (bounded), and for all references under
            # --all (dropdown / multi-plot / difference). exps drives the
            # cross-experiment selectors, so it is --all only.
            if single or embed_all:
                data_embed[slug] = _embed_data(grp, label, names, n_refs)
            if embed_all:
                exps.append({"slug": slug, "label": label, "n": n_refs})

        n_experiments = len(groups)

    if _can_compare(embed_all, exps):
        sections.extend(_comparison_sections())
    if structure:
        sections.append(
            _section(
                "Structure visualization", f'<div class="structure"><pre>{structure}</pre></div>'
            )
        )

    # The sequence view is client-rendered wherever a profile is (single-ref on
    # load, multi-ref via the dropdown), so embed it whenever arrays are embedded.
    seq_embed = sequences if (has_seq and data_embed) else None
    embed_tags = _embed_tags(figs_embed, data_embed, exps, seq_embed)
    return _document("\n".join(sections), embed_tags, n_experiments)


def main() -> None:
    """CLI entry point for ``cmuts plot``: render an HTML report from a reactivity h5."""
    parser = argparse.ArgumentParser(
        prog="cmuts-plot", description="Generate an HTML report from a reactivity HDF5 file"
    )
    parser.add_argument("file", help="HDF5 file written by cmuts normalize")
    parser.add_argument(
        "--all",
        action="store_true",
        help="Embed per-sequence plots, the multi-plot and the difference plot (small runs only)",
    )
    parser.add_argument("--group", default=None, help="Report only this group (default: all)")
    parser.add_argument("-o", "--out", default="report.html", help="Output HTML file")
    args = parser.parse_args()

    if not os.path.exists(args.file):
        raise FileNotFoundError(f"File not found: {args.file}")

    html = build(args.file, embed_all=args.all, group=args.group)
    with open(args.out, "w") as fh:
        fh.write(html)
    print(f"Report written to {args.out}")


if __name__ == "__main__":
    main()
