# `cmuts plot`

## Purpose

Serving an interactive HTML report summarizing and plotting the results contained in one or more `cmuts` HDF5 files.

## Requires

- One or more HDF5 files written by `cmuts hmm`, `sub`, `div`, or `norm`

As well as a Python 3 installation with the `h5py`, `numpy`, and `plotly` packages.

## Usage

The command starts a local web server and prints a link to open in your browser.

```sh
cmuts plot \
    subtracted.h5 --label "Subtracted" \
    treated.h5 --label "Treated" \
    untreated.h5 --label "Untreated"
```

The server runs until interrupted with Ctrl-C.

Saving it from the browser provides a static copy that can be viewed and shared offline. Only the currently selected reference is embedded, meaning reference selection is disabled in the static copy.

## Reactivity

The report computes the reactivity as the sum of the channels selected at the top of the page.

## Options

### Arguments

| Argument | Description |
| --- | --- |
| `HDF5` | outputs to show, one condition each |

### Input and output

| Option | Description |
| --- | --- |
| `--label NAME` | label for the file in the same position (default: the file's stem) |

### Serving

| Option | Description |
| --- | --- |
| `--host HOST` | host to bind (default 127.0.0.1) |
| `--port N` | port to bind (default: an OS-chosen free port) |

### Information

| Option | Description |
| --- | --- |
| `-h, --help` | show this help and exit |
| `-V, --version` | show the version and exit |
