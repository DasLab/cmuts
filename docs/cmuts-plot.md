# `cmuts plot`

## Purpose

Serving an interactive HTML report summarizing and plotting the results of one or more `cmuts` outputs.

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

The server runs until interrupted with Ctrl-C. The page embeds every condition's figures for the view it shows, so saving it from the browser gives a snapshot that can be viewed and shared offline: the condition chips still switch, while the reference dropdown and the channel checkboxes need the running server.

## Reactivity

The outputs hold one rate per event channel rather than one reactivity, so the report computes its reactivity as `cmuts norm` pools its scale: the aggregate of the mismatch, insertion, and deletion rates, one less the product of their no-event rates. Checkboxes at the top of the page select the channels the aggregate is taken over. Plotting an output of `cmuts norm` shows normalized profiles, since the norm divides every channel by the one scale.

## Conditions

The inputs need not hold the same library. Comparisons are by reference index and position, so a one-off construct can sit beside a large library, or a mutant beside its wild type. Each condition reports its own reference count, the sequence strip shows the selected condition's own bases, and a reference a condition does not hold simply shows nothing for it.

The report uses the combined single-reference layout only when every input holds one reference. Otherwise it shows aggregate and per-sequence sections, with the reference dropdown reaching the largest library.

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
