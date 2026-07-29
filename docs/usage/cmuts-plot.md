
## Purpose

`cmuts plot` renders a self-contained HTML report from a reactivity HDF5 file
written by `cmuts normalize`. The report contains, for every group in the file,
summary statistics and the aggregate plots (reactivity profiles, mutation
heatmap, terminations, coverage, read distributions, and the SNR-vs-read-depth
curves). Plotting is decoupled from normalization: normalize writes the data,
`cmuts plot` draws it, reading each dataset by name.

## Usage

```bash
cmuts plot INPUT.h5 -o report.html
```

Options:

- `-o`, `--out` — output HTML file (default `report.html`).
- `--group GROUP` — report only this group (default: every top-level group).
- `--all` — additionally embed every reference's profile and per-reference
  matrices (mutual information, correlation) behind a dropdown. This embeds
  per-reference data, so use it only for runs with one or a few references.

The report is a single HTML file with the plotting library inlined, so it opens
offline in any browser. To save a figure as an image, use its toolbar camera
button ("Download plot as PNG") — no extra tools are required.

An example can be found under `./examples/plot`.
