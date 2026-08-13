#!/usr/bin/env python3
"""Documentation generation based on the current binaries.

Programs describe their inputs and outputs via `--dump-layout` and
`--dump-options`. This script formats and inserts these data into the
documentation into sections marked by HTML markers.

    scripts/document.py build/release/cmuts-hmm docs/basics.md docs/cmuts-hmm.md
"""

import json
import subprocess
import sys
from pathlib import Path

# How a row is written, given the shape the program names and whether there is
# one row per reference. n is the references and l the longest of them; a field
# indexed by read length reaches twice that, a read being longer than what it
# aligns to.
EXTENTS = {"per base": "l", "per length": "2l", "scalar": None}


def table(headings: list, rows: list) -> str:
    """A markdown table. A cell holding a pipe would end its column early, so
    every cell is written with its pipes escaped."""
    def line(cells):
        return "| " + " | ".join(cell.replace("|", "\\|") for cell in cells) + " |"

    return "\n".join([line(headings), line("---" for _ in headings)]
                     + [line(row) for row in rows])


def described(program: str, flag: str) -> dict:
    """What a program says about itself.

    A program that will not answer is a page asking the wrong one, so the
    refusal is reported as it was given rather than raised through.
    """
    told = subprocess.run([program, flag], capture_output=True, text=True)

    if told.returncode != 0:
        raise SystemExit(f"{program} {flag}: {told.stderr.strip() or 'refused'}")

    return json.loads(told.stdout)


# ---------------------------------------------------------------------------
# The datasets an output holds
# ---------------------------------------------------------------------------


def shape(field: dict) -> str:
    """The dimensions of a field's dataset, as the documentation writes them."""
    if field["row"] not in EXTENTS:
        raise SystemExit(f"{field['name']}: unknown row shape {field['row']!r}")

    extents = [extent for extent in ("n" if field["per_reference"] else None,
                                     EXTENTS[field["row"]]) if extent]

    return "()" if not extents else f"({', '.join(extents)}{',' if len(extents) == 1 else ''})"


def layout(program: str) -> str:
    """The datasets an output holds, one to a row.

    The last column is the dataset's HDF5 fill value, which a reader can ask a
    file for. What seeing it means is the output page's, differing by field.
    """
    absent = {"nan": "NaN", "zero": "zero"}

    return table(
        ["Dataset", "Shape", "Type", "Fill"],
        [[f"`{field['name']}`", f"`{shape(field)}`", field["type"],
          absent.get(field["absent"], field["absent"])]
         for field in described(program, "--dump-layout")["fields"]],
    )


def fields(program: str) -> str:
    """The same datasets at length, each with what its numbers are.

    A field carrying no sentence is written without one rather than skipped, so
    that a dataset never goes unlisted for want of prose.
    """
    written = []

    for field in described(program, "--dump-layout")["fields"]:
        absent = "NaN" if field["absent"] == "nan" else field["absent"]
        written.append(f"`{field['name']}` — `{shape(field)}`, {field['type']}, "
                       f"fill {absent}.")

        if field["detail"]:
            written.append(field["detail"])

        written.append("")

    return "\n".join(written).rstrip()


# ---------------------------------------------------------------------------
# What the subtraction does to each of them
# ---------------------------------------------------------------------------


def rules(program: str) -> str:
    """What each dataset is combined by, with a denatured control and without.

    A field combined the same way either way is written once, since a column
    repeating its neighbour says only that the control changes nothing there.
    """
    rows = []

    for field in described(program, "--dump-rules")["fields"]:
        alone, controlled = field["uncontrolled"], field["controlled"]
        with_one = "the same" if controlled == alone else controlled["detail"]

        rows.append([f"`{field['name']}`", alone["detail"], with_one])

    return table(["Dataset", "Without a control", "With one"], rows)


# ---------------------------------------------------------------------------
# The arguments a program takes
# ---------------------------------------------------------------------------


def invocation(option: dict) -> str:
    """Both forms of an option and its placeholder, as the help prints them."""
    forms = [f"-{option['short']}"] if option["short"] else []
    forms.append(f"--{option['name']}")
    written = ", ".join(forms)

    return f"{written} {option['metavar']}" if option["metavar"] else written


def note(option: dict) -> str:
    """What the help says about an option beyond what it is for: the values it
    takes, and the one it has when it is not given.

    A range is given only where both ends are declared. Every count is bounded
    below by zero, which says nothing; a mapping quality stopping at 254 and a
    weight at 1 are the bounds worth reading.
    """
    notes = []

    if option["choices"]:
        notes.append(", ".join(option["choices"]))

    if option["minimum"] is not None and option["maximum"] is not None:
        notes.append(f"{option['minimum']} to {option['maximum']}")

    if option["required"]:
        notes.append("required")
    elif option["unset_label"]:
        notes.append(f"default: {option['unset_label']}")
    elif option["type"] != "flag" and option["default"] is not None:
        notes.append(f"default {option['default']}")

    return f" ({'; '.join(notes)})" if notes else ""


def options(program: str) -> str:
    """Every argument a program takes, under the headings its help groups them
    by. The hidden ones describe the program to a machine and are left out."""
    spoken = described(program, "--dump-options")
    written = []

    if spoken["positionals"]:
        rows = [[f"`{p['metavar']}{'...' if p['variadic'] else ''}`", p["help"]]
                for p in spoken["positionals"]]
        written += ["### Arguments", "", table(["Argument", "Description"], rows), ""]

    shown = [option for option in spoken["options"] if not option["hidden"]]

    for group in dict.fromkeys(option["group"] for option in shown):
        rows = [[f"`{invocation(option)}`", option["help"] + note(option)]
                for option in shown if option["group"] == group]
        written += [f"### {group}", "", table(["Option", "Description"], rows), ""]

    return "\n".join(written).rstrip()


# ---------------------------------------------------------------------------
# Writing them into the pages
# ---------------------------------------------------------------------------

WRITERS = {"LAYOUT": layout, "FIELDS": fields, "RULES": rules, "OPTIONS": options}


def spliced(text: str, kind: str, generated: str) -> str:
    """The page with what lies between one pair of markers replaced."""
    begin, end = f"<!-- BEGIN GENERATED {kind} -->", f"<!-- END GENERATED {kind} -->"
    start, stop = text.find(begin), text.find(end)

    if start < 0:
        return text

    if stop < 0:
        raise SystemExit(f"{begin} with no {end} after it")

    return text[:start] + begin + "\n" + generated + "\n" + text[stop:]


def main(program: str, pages: list) -> None:
    for name in pages:
        page = Path(name)
        text = page.read_text()
        written = [kind for kind in WRITERS if f"<!-- BEGIN GENERATED {kind} -->" in text]

        if not written:
            raise SystemExit(f"{page}: no generated block to write into")

        for kind in written:
            text = spliced(text, kind, WRITERS[kind](program))

        page.write_text(text)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)

    main(sys.argv[1], sys.argv[2:])
