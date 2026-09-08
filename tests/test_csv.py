"""Writing an output as comma separated values.

The table depends only on the datasets of the input, so the files it is run over
are written by hand. The sequence supplies the base at each position and the
length of each reference. A FASTA, where one is given, supplies the names only.
"""

from __future__ import annotations

import csv
import io
import itertools

import numpy as np
import pytest

from alignments import write_fasta
from inputs import CAP, N_REFS, not_hdf5, random_fields
from oracle import sequences
from outputs import (
    COVERAGE,
    DELETION_ERROR,
    DELETION_RATE,
    INSERTION_ERROR,
    INSERTION_RATE,
    MISMATCH_ERROR,
    MISMATCH_RATE,
    SEQUENCE,
    delete_field,
    field_of,
)
from programs import CMUTS_CSV, attempt, run_cmuts, run_csv, try_csv

# The columns identifying a row, and the datasets that follow them, in the order
# they are written.
KEY_COLUMNS = ("reference", "position", "base")
VALUE_COLUMNS = (MISMATCH_RATE, MISMATCH_ERROR, INSERTION_RATE, INSERTION_ERROR,
                 DELETION_RATE, DELETION_ERROR, COVERAGE)

# The character of each token, and the value a column past the end of a
# reference holds.
TOKENS = "ACGTN"
OUTSIDE = -1

# The bases the build fixture writes into every row, which cycle through the four
# named bases across a row of CAP columns.
BASES = "".join(TOKENS[at % 4] for at in range(CAP))

# One name for each row of a file the build fixture writes.
NAMES = ("a", "b", "c", "d")

# Values are written to six significant digits, so a comparison against the
# dataset they came from allows for the rounding.
TOLERANCE = 1e-5


@pytest.fixture
def named(tmp_path):
    """Returns a function that writes an output as a table, naming the rows from
    a FASTA written for it. Each FASTA is named separately, so one test may
    convert several."""
    written = itertools.count()

    def run(path, references, **options):
        fasta = write_fasta(references, tmp_path / f"refs{next(written)}.fasta")

        return run_csv(path, fasta, **options)

    return run


def every_reference(bases: str = BASES) -> dict:
    """One record per row of a file the build fixture writes."""
    return {name: bases for name in NAMES}


def rows_holding(*given) -> np.ndarray:
    """Builds a sequence dataset from the bases of each reference in turn, each
    padded out to the width of a row. A row the caller does not give holds
    padding alone, and so holds no reference."""
    rows = np.full((N_REFS, CAP), OUTSIDE, dtype=np.int8)

    for row, bases in enumerate(given):
        rows[row, :len(bases)] = [TOKENS.index(base) for base in bases]

    return rows


def table_of(text: str) -> tuple:
    """Returns the header and the rows of a table, each row as its cells. The
    reader takes quoted fields apart, so a test reads a name as it was
    written."""
    lines = list(csv.reader(io.StringIO(text)))

    assert lines

    return lines[0], lines[1:]


def column_of(header, rows, name, reference) -> list:
    """Returns one column of one reference's rows."""
    at = header.index(name)

    return [row[at] for row in rows if row[0] == reference]


# ---------------------------------------------------------------------------
# The columns
# ---------------------------------------------------------------------------


def test_the_header_names_the_key_columns_and_every_dataset(build, named):
    header, _ = table_of(named(build(), every_reference()))

    assert header == [*KEY_COLUMNS, *VALUE_COLUMNS]


@pytest.mark.parametrize("missing", [INSERTION_RATE, DELETION_ERROR])
def test_a_dataset_the_input_lacks_gets_no_column(build, named, missing):
    header, _ = table_of(named(delete_field(build(), missing), every_reference()))

    assert missing not in header


# ---------------------------------------------------------------------------
# The rows
# ---------------------------------------------------------------------------


def test_every_position_of_every_reference_gets_a_row(build, named):
    _, rows = table_of(named(build(), every_reference()))

    assert len(rows) == N_REFS * CAP


def test_the_rows_follow_the_file(build, named):
    _, rows = table_of(named(build(), every_reference()))

    assert [row[0] for row in rows][::CAP] == list(NAMES)


def test_the_positions_count_from_one(build, named):
    header, rows = table_of(named(build(), every_reference()))

    assert column_of(header, rows, "position", NAMES[0]) == \
        [str(at + 1) for at in range(CAP)]


# ---------------------------------------------------------------------------
# The sequence
# ---------------------------------------------------------------------------


def test_the_bases_come_from_the_input(build, named):
    """The FASTA supplies the names only, so the base column is the input's own
    sequence even where the two differ."""
    other = "GGGGGG"
    header, rows = table_of(named(build(), every_reference(other)))

    assert other != BASES
    assert "".join(column_of(header, rows, "base", NAMES[0])) == BASES


@pytest.mark.parametrize("bases", ["ACGT", "ACGNN", "NNNNNN"])
def test_a_row_ends_where_its_reference_does(build, bases):
    """A column past the end of a reference holds padding and gets no row. The
    token for N is not that padding, so a reference ending in N keeps its full
    length."""
    header, rows = table_of(run_csv(build({SEQUENCE: rows_holding(bases)})))

    assert len(rows) == len(bases)
    assert "".join(column_of(header, rows, "base", "1")) == bases


def test_a_row_holding_no_reference_gets_no_row(build):
    _, rows = table_of(run_csv(build({SEQUENCE: rows_holding("", "AC")})))

    assert len(rows) == len("AC")
    assert [row[0] for row in rows] == ["2"] * len("AC")


# ---------------------------------------------------------------------------
# The names
# ---------------------------------------------------------------------------


def test_the_rows_are_numbered_from_one_without_a_fasta(build):
    _, rows = table_of(run_csv(build()))

    assert [row[0] for row in rows][::CAP] == [str(row + 1) for row in range(N_REFS)]


def test_a_name_holding_a_separator_is_quoted(build, named):
    """An unquoted comma or quote inside a name would split it across cells, so
    the name must come back whole."""
    name = 'a,b"c'

    _, rows = table_of(named(build(n_refs=1), {name: BASES}))

    assert all(row[0] == name for row in rows)


# ---------------------------------------------------------------------------
# The values
# ---------------------------------------------------------------------------


def test_each_value_is_the_dataset_at_that_position(build, named):
    values = random_fields(seed=61)
    values[SEQUENCE] = rows_holding(*([BASES] * N_REFS))
    header, rows = table_of(named(build(values), every_reference()))

    for name in VALUE_COLUMNS:
        written = [float(cell) for cell in column_of(header, rows, name, NAMES[1])]

        assert np.allclose(written, values[name][1], rtol=TOLERANCE), name


def test_a_value_that_is_not_a_number_is_an_empty_field(build, named):
    header, rows = table_of(named(build({MISMATCH_RATE: np.nan}), every_reference()))

    assert column_of(header, rows, MISMATCH_RATE, NAMES[0]) == [""] * CAP


# ---------------------------------------------------------------------------
# What is refused
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("count", [N_REFS - 1, N_REFS + 1])
def test_a_fasta_naming_another_number_of_references_is_refused(build, tmp_path, count):
    records = {f"ref{at}": BASES for at in range(count)}
    fasta = write_fasta(records, tmp_path / "refs.fasta")

    assert try_csv(build(), fasta).returncode == 1


def test_a_record_of_another_length_is_refused(build, tmp_path):
    """A record of another length describes another molecule, so the FASTA is not
    the one the input was counted against."""
    fasta = write_fasta(every_reference(BASES[:-1]), tmp_path / "refs.fasta")

    failed = try_csv(build(), fasta)

    assert failed.returncode == 1
    assert "not the one it was counted against" in failed.stderr


@pytest.mark.parametrize("row", [
    [0, 9, 2, OUTSIDE, OUTSIDE, OUTSIDE],
    [0, -3, 2, OUTSIDE, OUTSIDE, OUTSIDE],
    [0, 1, OUTSIDE, 2, OUTSIDE, OUTSIDE],
])
def test_a_sequence_that_is_not_tokens_then_padding_is_refused(build, row):
    """No run writes such a row. Left unchecked, a value that is not a token
    would print as a NUL byte, and a token after the padding would be dropped
    along with every value after it."""
    rows = np.tile(np.array(row, dtype=np.int8), (N_REFS, 1))

    failed = try_csv(build({SEQUENCE: rows}))

    assert failed.returncode == 1
    assert "not a sequence" in failed.stderr


@pytest.mark.parametrize("missing", [COVERAGE, SEQUENCE])
def test_an_input_without_a_required_dataset_is_refused(build, missing):
    assert try_csv(delete_field(build(), missing)).returncode == 1


def test_an_input_that_is_not_an_output_is_refused(tmp_path):
    assert try_csv(not_hdf5(tmp_path)).returncode == 1


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------


def test_the_input_is_required():
    assert attempt([*CMUTS_CSV]).returncode == 2


# ---------------------------------------------------------------------------
# End to end
# ---------------------------------------------------------------------------


def test_cmuts_csv_reads_what_cmuts_hmm_writes(data, falsifiable, tmp_path):
    rates = run_cmuts(data, tmp_path / "rates.h5")
    references = sequences(data.fasta)

    falsifiable(len(references) > 0)

    header, rows = table_of(run_csv(rates, data.fasta))

    assert header[:len(KEY_COLUMNS)] == list(KEY_COLUMNS)
    assert len(rows) == sum(len(bases) for bases in references.values())
    assert "".join(row[2] for row in rows) == "".join(references.values()).upper()


def test_the_table_holds_the_sequence_the_input_holds(data, falsifiable, tmp_path):
    """The base column is decoded from the sequence dataset, so it matches that
    dataset read directly."""
    rates = run_cmuts(data, tmp_path / "rates.h5")
    written = field_of(rates, SEQUENCE)

    falsifiable(written.shape[0] > 0)

    _, rows = table_of(run_csv(rates))
    tokens = [row[row >= 0] for row in written]

    assert "".join(row[2] for row in rows) == \
        "".join(TOKENS[token] for row in tokens for token in row)
