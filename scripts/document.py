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

# The extents a row occupies, given the shape the program names. n is the
# references and l the longest of them; a field indexed by read length reaches
# twice that, a read being longer than what it aligns to, and one indexed by a
# pair of positions is square.
EXTENTS = {"per base": ("l",), "per length": ("2l",), "per pair": ("l", "l"),
           "scalar": ()}

# Whether a run writes a field at all, for the column that says so.
PRESENT = {False: "always", True: "when asked for"}

# What marks a field a run may leave out of its output.
OPTIONAL = "written only when asked for"

# The value a dataset holds where the run wrote nothing, spelled as the pages
# spell it. What it means differs by field, and the output page covers it.
ABSENT = {"nan": "NaN", "zero": "zero"}

# The class each field's heading carries, which the stylesheet draws a rule
# beside. MyST attaches it to the section the heading opens, so it needs the
# attrs_block extension that docs/conf.py enables.
FIELD_CLASS = "{.field}"


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

    extents = ["n"] * field["per_reference"] + list(EXTENTS[field["row"]])

    if not extents:
        return "()"

    return f"({', '.join(extents)}{',' if len(extents) == 1 else ''})"


def layout(program: str) -> str:
    """The datasets an output holds, one to a row.

    The last column is the dataset's HDF5 fill value, which h5py reports as
    dataset.fillvalue.
    """
    return table(
        ["Dataset", "Shape", "Type", "Fill", "Written"],
        [[f"`{field['name']}`", f"`{shape(field)}`", f"`{field['type']}`",
          f"`{ABSENT.get(field['absent'], field['absent'])}`",
          PRESENT[field["optional"]]]
         for field in described(program, "--dump-layout")["fields"]],
    )


def attributes(program: str) -> str:
    """The attributes an output carries, one to a row.

    These sit on the root group and describe the run, so h5py reads them from
    the file's own attrs and not from any dataset.
    """
    return table(
        ["Attribute", "Description"],
        [[f"`{attribute['name']}`", attribute["detail"]]
         for attribute in described(program, "--dump-layout")["attributes"]],
    )


def fields(program: str) -> str:
    """The same datasets at length, each under a heading of its own.

    The heading gives every dataset a permalink and an entry in the table of
    contents, and leaves the description as much room as it needs. A field
    with no description gets a heading anyway, so that the block covers every
    dataset before the prose is written.

    Each heading opens a section of its own in the rendered page and carries
    the class the stylesheet draws a rule beside.
    """
    written = []

    for field in described(program, "--dump-layout")["fields"]:
        absent = ABSENT.get(field["absent"], field["absent"])
        marked = f" · _{OPTIONAL}_" if field["optional"] else ""
        written += [FIELD_CLASS, f"### `{field['name']}`", "",
                    f"**Shape** `{shape(field)}` · **Type** `{field['type']}` · "
                    f"**Fill** `{absent}`{marked}", ""]

        if field["detail"]:
            written += [field["detail"], ""]

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


def option_rows(options: list) -> list:
    """The rows a table of options holds."""
    return [[f"`{invocation(option)}`", option["help"] + note(option)]
            for option in options]


def options(program: str) -> str:
    """Every argument a program takes, under the headings the help groups them
    by.

    An option the help leaves out goes under Advanced, so that a page describes
    everything the program accepts and the help stays short.
    """
    spoken = described(program, "--dump-options")
    written = []

    if spoken["positionals"]:
        rows = [[f"`{p['metavar']}{'...' if p['variadic'] else ''}`", p["help"]]
                for p in spoken["positionals"]]
        written += ["### Arguments", "", table(["Argument", "Description"], rows), ""]

    shown  = [option for option in spoken["options"] if not option["hidden"]]
    hidden = [option for option in spoken["options"] if option["hidden"]]

    for group in dict.fromkeys(option["group"] for option in shown):
        rows = option_rows([option for option in shown if option["group"] == group])
        written += [f"### {group}", "", table(["Option", "Description"], rows), ""]

    if hidden:
        written += ["### Advanced", "",
                    "Accepted, and left out of `--help`.", "",
                    table(["Option", "Description"], option_rows(hidden)), ""]

    return "\n".join(written).rstrip()


# ---------------------------------------------------------------------------
# Writing them into the pages
# ---------------------------------------------------------------------------

WRITERS = {"LAYOUT": layout, "ATTRIBUTES": attributes, "FIELDS": fields,
           "OPTIONS": options}

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
