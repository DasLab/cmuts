#!/usr/bin/env python3
"""Documentation generation based on the current binaries.

Programs describe their inputs and outputs via `--dump-layout` and
`--dump-options`. This script formats and inserts these data into the
documentation into sections marked by HTML markers.

    scripts/document.py build/release docs
"""

import json
import re
import subprocess
import sys
from pathlib import Path

# How a row is written, given the shape the program names and whether there is
# one row per reference. n is the references and l the longest of them; a field
# indexed by read length reaches twice that, a read being longer than what it
# aligns to.
EXTENTS = {"per base": "l", "per length": "2l", "scalar": None}

# The value a dataset holds where the run wrote nothing, spelled as the pages
# spell it. What it means differs by field, and the output page covers it.
ABSENT = {"nan": "NaN", "zero": "zero"}


def table(headings: list, rows: list) -> str:
    """A markdown table. A cell holding a pipe would end its column early, so
    every cell is written with its pipes escaped."""
    def line(cells):
        return "| " + " | ".join(cell.replace("|", "\\|") for cell in cells) + " |"

    return "\n".join([line(headings), line("---" for _ in headings)]
                     + [line(row) for row in rows])


def described(program: str, flag: str) -> dict:
    """Runs one of a program's dump flags and parses the JSON it prints.

    A non-zero exit usually means the page asked for a dump this program does
    not have, so the message it printed is shown instead of a traceback.
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

    The last column is the dataset's HDF5 fill value, which h5py reports as
    dataset.fillvalue.
    """
    return table(
        ["Dataset", "Shape", "Type", "Fill"],
        [[f"`{cell}`" for cell in (field["name"], shape(field), field["type"],
                                   ABSENT.get(field["absent"], field["absent"]))]
         for field in described(program, "--dump-layout")["fields"]],
    )


def fields(program: str) -> str:
    """The same datasets at length, each under a heading of its own.

    The heading gives every dataset a permalink and an entry in the table of
    contents, and leaves the description as much room as it needs. A field
    with no description gets a heading anyway, so that the block covers every
    dataset before the prose is written.

    Each one is wrapped in a div the stylesheet draws a rule beside. The blank
    lines around the tags are what keep the markdown between them markdown,
    both here and on GitHub, which has no such stylesheet and shows the fields
    plainly.
    """
    written = []

    for field in described(program, "--dump-layout")["fields"]:
        absent = ABSENT.get(field["absent"], field["absent"])
        written += ['<div class="field" markdown>', "",
                    f"### `{field['name']}`", "",
                    f"**Shape** `{shape(field)}` · **Type** `{field['type']}` · "
                    f"**Fill** `{absent}`", ""]

        if field["detail"]:
            written += [field["detail"], ""]

        written += ["</div>", ""]

    return "\n".join(written).rstrip()


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
    """The values an option takes and the one it has when it is not given.

    A range is included only where both ends are declared. Every count has a
    floor of zero, which is not worth printing; the useful ranges are the
    mapping quality's 0 to 254 and a weight's 0 to 1.
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
    """Every argument a program takes, under the headings the help groups them
    by. Hidden options exist for this script to read and are left out."""
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

WRITERS = {"LAYOUT": layout, "FIELDS": fields, "OPTIONS": options}

# A block names the program that answers it and the table wanted from it, so a
# page is found by looking in it rather than by being listed somewhere.
BLOCK = re.compile(r"<!-- BEGIN GENERATED (?P<program>[\w.-]+) (?P<kind>[A-Z]+) -->")


def spliced(text: str, program: str, kind: str, generated: str) -> str:
    """The page with what lies between one pair of markers replaced."""
    begin = f"<!-- BEGIN GENERATED {program} {kind} -->"
    end = f"<!-- END GENERATED {program} {kind} -->"
    start, stop = text.find(begin), text.find(end)

    if stop < 0:
        raise SystemExit(f"{begin} with no {end} after it")

    return text[:start] + begin + "\n" + generated + "\n" + text[stop:]


def main(build: str, docs: str) -> None:
    for page in sorted(Path(docs).glob("*.md")):
        text = page.read_text()
        blocks = BLOCK.findall(text)

        for program, kind in blocks:
            binary = Path(build) / program

            if kind not in WRITERS:
                raise SystemExit(f"{page}: no table named {kind}")

            if not binary.exists():
                raise SystemExit(f"{page}: {binary} is not built")

            text = spliced(text, program, kind, WRITERS[kind](str(binary)))

        if blocks:
            page.write_text(text)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)

    main(sys.argv[1], sys.argv[2])
