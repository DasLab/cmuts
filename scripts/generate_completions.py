#!/usr/bin/env python3
"""Generate bash and zsh completion scripts for cmuts using shtab.

cmuts is a bash dispatcher over subcommands implemented in different languages:
C++ (``core``, ``generate``; arguments declared in ``config/*.json``), bash
(``align``, ``test``), and Python argparse (``normalize``, ``plot``,
``visualize``). This module assembles a single argparse parser describing the
whole CLI — importing the Python subcommands' parsers live so their completions
can never drift from the argparse definitions — and hands it to shtab to emit
the shell scripts.
"""

from __future__ import annotations

import argparse
import importlib.machinery
import importlib.util
import json
import sys
from pathlib import Path

try:
    import shtab
except ImportError:  # pragma: no cover - completions are optional at build time
    shtab = None  # type: ignore[assignment]

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / "config"
SRC_PY = ROOT / "src" / "python"
COMPLETIONS = ROOT / "completions"


def _file(action: argparse.Action, *exts: str) -> argparse.Action:
    """Complete file paths for ``action`` (filtered by extension under zsh)."""
    if exts:
        glob = "|".join(exts)
        action.complete = {"bash": "_shtab_compgen_files", "zsh": f"_files -g '*.({glob})'"}
    else:
        action.complete = shtab.FILE
    return action


def _annotate_files(parser: argparse.ArgumentParser, mapping: dict[str, tuple[str, ...]]) -> None:
    """Set file completion on each action matching an option string or dest."""
    for action in parser._actions:
        keys = set(action.option_strings) | {action.dest}
        for key, exts in mapping.items():
            if key in keys:
                _file(action, *exts)


def _load_script_parser(script: str) -> argparse.ArgumentParser:
    """Import the module-level ``parser`` from an extensionless CLI script."""
    name = "_cmuts_cli_" + script.replace("-", "_")
    loader = importlib.machinery.SourceFileLoader(name, str(SRC_PY / script))
    spec = importlib.util.spec_from_loader(name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module.parser


def _json_subparser(sub, name: str, json_file: str, help_text: str) -> argparse.ArgumentParser:
    """Build a subparser for a C++ program from its argument JSON config."""
    with open(CONFIG / json_file) as f:
        config = json.load(f)
    p = sub.add_parser(name, help=help_text)
    for group in config["groups"]:
        for arg in group["args"]:
            if not arg["long"].startswith("--"):
                continue
            names = ([arg["short"]] if arg.get("short") else []) + [arg["long"]]
            kwargs: dict = {"help": arg.get("help", "")}
            if arg["type"] == "bool":
                kwargs["action"] = "store_true"
            p.add_argument(*names, **kwargs)
    return p


def build_parser() -> argparse.ArgumentParser:
    """Assemble the full cmuts CLI as one argparse parser for shtab."""
    parser = argparse.ArgumentParser(prog="cmuts")
    parser.add_argument("-v", "--version", action="store_true", help="Show version")
    _file(parser.add_argument("--log", help="Log file"))
    sub = parser.add_subparsers(dest="subcommand")

    # --- C++ subcommands: arguments come from config/*.json ---
    core = _json_subparser(
        sub, "core", "cmuts_args.json", "Count chemical mutations from aligned reads"
    )
    core.add_argument("--threads", help="Number of MPI processes")
    _file(core.add_argument("inputs", nargs="*", help="Input alignments"), "sam", "bam", "cram")
    _annotate_files(core, {"--fasta": ("fa", "fasta", "fna"), "--output": ("h5",)})

    gen = _json_subparser(sub, "generate", "tests_args.json", "Generate synthetic test data")
    _annotate_files(
        gen,
        {
            "--out-fasta": ("fa", "fasta", "fna"),
            "--out": ("sam", "bam", "cram"),
            "--out-h5": ("h5",),
        },
    )

    # --- Python subcommands: parsers imported live so completions never drift ---
    from cmuts.normalize.__main__ import parser as normalize_parser

    norm = sub.add_parser(
        "normalize",
        help="Normalize mutation counts to reactivity profiles",
        add_help=False,
        parents=[normalize_parser],
    )
    _annotate_files(norm, {"file": ("h5",), "--fasta": ("fa", "fasta", "fna"), "--out": ("h5",)})

    # cmuts-plot delegates to report:main (no module-level parser to import).
    plot = sub.add_parser("plot", help="Render an HTML report from a normalize h5")
    _file(plot.add_argument("file", help="Reactivity HDF5 file"), "h5")
    plot.add_argument(
        "--all", action="store_true", help="Embed per-reference plots behind a dropdown"
    )
    plot.add_argument("--group", help="Report only this group")
    _file(plot.add_argument("-o", "--out", help="Output HTML file"), "html")

    vis = sub.add_parser(
        "visualize",
        help="Color 3D structures by reactivity",
        add_help=False,
        parents=[_load_script_parser("cmuts-visualize")],
    )
    _annotate_files(
        vis, {"--file": ("h5",), "--cif": ("cif", "mmcif"), "--fasta": ("fa", "fasta", "fna")}
    )

    # --- bash subcommands: declared here to match their dispatchers ---
    align = sub.add_parser("align", help="Align FASTQ files to a reference")
    _file(align.add_argument("--fasta", help="Reference FASTA"), "fa", "fasta", "fna")
    _file(align.add_argument("--output", help="Output directory"))
    align.add_argument("--overwrite", action="store_true", help="Overwrite existing output")
    align.add_argument("--rebuild", action="store_true", help="Rebuild index files")
    align.add_argument("--verbose", action="store_true", help="Verbose output")
    align.add_argument("--threads", help="Number of threads")
    _file(align.add_argument("--barcodes", help="Barcode file"), "csv", "tsv", "txt")
    align.add_argument("--trim-5", help="5' adapter sequence")
    align.add_argument("--trim-3", help="3' adapter sequence")
    align.add_argument("--local", action="store_true", help="Local alignment mode")
    _file(align.add_argument("--pairs", help="Paired-end reads"), "fastq", "fastq.gz", "fq", "fq.gz")
    _file(
        align.add_argument("inputs", nargs="*", help="Input reads"),
        "fastq",
        "fastq.gz",
        "fq",
        "fq.gz",
        "bam",
    )

    test = sub.add_parser("test", help="Run integration tests")
    test.add_argument("--quick", action="store_true", help="Run quick tests only")
    test.add_argument("--cram", action="store_true", help="Run CRAM-specific tests only")

    return parser


def main() -> None:
    if shtab is None:
        print("shtab not installed; skipping completion regeneration.", file=sys.stderr)
        return
    try:
        parser = build_parser()
    except ImportError as exc:
        print(f"cannot import cmuts CLI parsers ({exc}); skipping completions.", file=sys.stderr)
        return
    COMPLETIONS.mkdir(exist_ok=True)
    (COMPLETIONS / "cmuts.bash").write_text(shtab.complete(parser, "bash"))
    (COMPLETIONS / "_cmuts").write_text(shtab.complete(parser, "zsh"))
    print(f"Generated {COMPLETIONS / 'cmuts.bash'}")
    print(f"Generated {COMPLETIONS / '_cmuts'}")


if __name__ == "__main__":
    main()
