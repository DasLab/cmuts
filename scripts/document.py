#!/usr/bin/env python3
"""Writes the parts of the documentation that describe the programs.

A program describes itself, so nothing here restates what one holds: the table
of datasets comes from `cmuts-hmm --dump-layout`, which reads the declaration
the program writes its output from.

Only what lies between a pair of markers is replaced, so a page is edited
around them and never inside them.

    scripts/document.py build/release/cmuts-hmm docs/basics.md
"""

import json
import subprocess
import sys
from pathlib import Path

BEGIN = "<!-- BEGIN GENERATED LAYOUT -->"
END = "<!-- END GENERATED LAYOUT -->"

# How a row is written, given the shape the program names and whether there is
# one row per reference. n is the references and l the longest of them; a field
# indexed by read length reaches twice that, a read being longer than what it
# aligns to.
EXTENTS = {"per base": "l", "per length": "2l", "scalar": None}


def shape(field: dict) -> str:
    """The dimensions of a field's dataset, as the documentation writes them."""
    if field["row"] not in EXTENTS:
        raise SystemExit(f"{field['name']}: unknown row shape {field['row']!r}")

    extents = [extent for extent in ("n" if field["per_reference"] else None,
                                     EXTENTS[field["row"]]) if extent]

    return "()" if not extents else f"({', '.join(extents)}{',' if len(extents) == 1 else ''})"


def table(layout: dict) -> str:
    """The datasets an output holds, one to a row.

    The last column is the dataset's HDF5 fill value, which a reader can ask a
    file for. What seeing it means is the output page's, differing by field.
    """
    header = ["| Dataset | Shape | Type | Fill |", "| --- | --- | --- | --- |"]
    absent = {"nan": "NaN", "zero": "zero"}

    rows = [
        f"| `{field['name']}` | `{shape(field)}` | {field['type']} "
        f"| {absent.get(field['absent'], field['absent'])} |"
        for field in layout["fields"]
    ]

    return "\n".join(header + rows)


def spliced(text: str, generated: str, page: Path) -> str:
    """The page with what lies between the markers replaced."""
    start, end = text.find(BEGIN), text.find(END)

    if start < 0 or end < 0:
        raise SystemExit(f"{page}: no {BEGIN} ... {END} to write into")

    return text[:start] + BEGIN + "\n" + generated + "\n" + text[end:]


def main(program: str, pages: list) -> None:
    described = subprocess.run([program, "--dump-layout"], check=True,
                               capture_output=True, text=True).stdout
    generated = table(json.loads(described))

    for name in pages:
        page = Path(name)
        page.write_text(spliced(page.read_text(), generated, page))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)

    main(sys.argv[1], sys.argv[2:])
