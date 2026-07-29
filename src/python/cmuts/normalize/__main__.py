#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import sys

import h5py

import cmuts
from cmuts.internal import _names_from_fasta

NAME = "cmuts normalize"


parser = argparse.ArgumentParser()
parser.add_argument(
    "file",
    help="The HDF5 file containing the input data.",
)
parser.add_argument(
    "--experiment",
    action="append",
    nargs="+",
    required=True,
    metavar="NAME mod=... [nomod=...]",
    help=(
        "Define an experiment: a NAME followed by 'mod=' and optional 'nomod=' "
        "(comma-separated HDF5 paths; replicates are summed). The NAME is the "
        "output HDF5 group. Repeat the flag for multiple experiments, e.g. "
        "--experiment apo mod=a.h5,b.h5 nomod=c.h5 --experiment holo mod=d.h5 nomod=e.h5"
    ),
)
parser.add_argument(
    "--per-experiment-norm",
    action="store_true",
    help=(
        "Normalize each experiment independently. By default one factor is "
        "shared across all experiments, keeping them on a comparable scale."
    ),
)
parser.add_argument(
    "--per-reference-norm",
    action="store_true",
    help=(
        "Normalize each reference independently (like rf-norm per transcript). "
        "By default one factor is shared across all references."
    ),
)
parser.add_argument(
    "--fasta",
    help="FASTA file of reference sequences.",
    required=True,
)
parser.add_argument(
    "-o", "--out", help="The output HDF5 file to write to.", default="reactivity.h5"
)
parser.add_argument(
    "--overwrite",
    help="Overwrite an existing HDF5 file (the whole file, not just the group).",
    action="store_true",
)
parser.add_argument(
    "--clip-below",
    help="Clip reactivity values below this threshold (e.g. 0).",
    type=float,
    default=None,
    metavar="N",
)
parser.add_argument(
    "--clip-above",
    help="Clip reactivity values above this threshold (e.g. 1).",
    type=float,
    default=None,
    metavar="N",
)
parser.add_argument(
    "--blank-5p",
    help="NaN out this many bases on the 5' end.",
    type=int,
    default=0,
)
parser.add_argument(
    "--blank-3p",
    help="NaN out this many bases on the 3' end.",
    type=int,
    default=0,
)
parser.add_argument(
    "--blank-cutoff",
    help="NaN out any positions with less than this many reads.",
    type=int,
    default=10,
)
parser.add_argument(
    "--norm",
    help="Normalization scheme (see the cmuts normalize docs for what each does).",
    choices=cmuts.normalize.scheme_names(),
    default="ubr",
)
parser.add_argument(
    "--no-insertions",
    help="Do not use insertions to compute the reactivity profile.",
    action="store_true",
)
parser.add_argument(
    "--no-deletions",
    help="Do not use deletions to compute the reactivity profile.",
    action="store_true",
)
parser.add_argument(
    "--sig",
    help="The significance level for correlation plots.",
    type=float,
    default=0.05,
)


def _parse_experiment(tokens: list[str]) -> cmuts.Experiment:
    """Parse one ``--experiment NAME mod=... [nomod=...]`` token list.

    The first token is the experiment name; the rest are ``key=value`` pairs
    (``mod`` required, ``nomod`` optional) whose values are comma-separated HDF5
    paths.
    """
    if not tokens:
        sys.exit("--experiment requires a NAME followed by mod=...")
    name, rest = tokens[0], tokens[1:]
    mod: list[str] | None = None
    nomod: list[str] | None = None
    for tok in rest:
        if "=" not in tok:
            sys.exit(f"--experiment {name}: expected key=value, got '{tok}'.")
        key, _, value = tok.partition("=")
        paths = [p for p in value.split(",") if p]
        if key == "mod":
            mod = paths
        elif key == "nomod":
            nomod = paths
        else:
            sys.exit(f"--experiment {name}: unknown key '{key}' (use mod= or nomod=).")
    if not mod:
        sys.exit(f"--experiment {name}: needs a non-empty mod=...")
    return cmuts.Experiment(name=name, mod=mod, nomod=nomod)


def _remove_if_exists(path: str, overwrite: bool = False) -> None:
    if os.path.exists(path):
        if overwrite:
            os.remove(path)


def main():
    args = parser.parse_args()

    experiments = [_parse_experiment(tokens) for tokens in args.experiment]
    names = [e.name for e in experiments]
    if len(set(names)) != len(names):
        parser.error("--experiment names must be unique (each is an output HDF5 group).")
    if "meta" in names:
        parser.error("'meta' is a reserved group name (used for shared metadata).")

    print()
    print(cmuts.title(NAME, cmuts.__version__))
    print(cmuts.subtitle(", ".join(names)))

    _remove_if_exists(args.out, args.overwrite)

    clip = (args.clip_below, args.clip_above)
    blank = (args.blank_5p, args.blank_3p)

    # opts.mod and opts.nomod are unused by compute_reactivities; pass empties.
    opts = cmuts.Opts(
        cmuts.DataGroups([]),
        cmuts.DataGroups(None),
        args.blank_cutoff,
        not args.no_insertions,
        not args.no_deletions,
        args.norm,
        blank,
        clip,
        args.sig,
        per_experiment=args.per_experiment_norm,
        per_reference=args.per_reference_norm,
    )

    with h5py.File(args.file, "r") as f:
        results = cmuts.compute_reactivities(f, args.fasta, experiments, opts)

    # Save all experiments to a single output file (one HDF5 group per name).
    # Reference names (FASTA order) are stored so the report can label by name.

    names = _names_from_fasta(args.fasta)
    cmuts.save_groups(args.out, [(r.experiment.name, r.combined) for r in results], names=names)

    # Print stats and save the SNR-vs-depth curves alongside each experiment.
    # Plotting is decoupled: `cmuts plot <out.h5>` renders the figures.
    for r in results:
        if len(results) > 1:
            print()
            print(cmuts.subtitle(r.experiment.name))
        cmuts.stats(r.mod, r.nomod, r.combined)
        cmuts.compute_snr_curves(r.mod, r.nomod, r.combined).save(r.experiment.name, args.out)


if __name__ == "__main__":
    main()
