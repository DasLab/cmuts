"""What happens to a file already at the output path.

Replacing a previous result requires --overwrite, a run costing far more than
the command that starts it.
"""

import pytest

from support import outputs_agree, run_cmuts, try_cmuts


def test_an_existing_output_is_not_replaced_without_overwrite(data, tmp_path):
    output = tmp_path / "out.h5"
    run_cmuts(data, output)
    before = output.read_bytes()

    attempt = try_cmuts(data, output)

    assert attempt.returncode != 0
    assert "already holds data" in attempt.stderr
    assert output.read_bytes() == before, "the first result is untouched"


def test_overwrite_replaces_an_existing_output(data, tmp_path):
    """The two runs are given different criteria, so that agreeing with one
    identifies which result was left at the path."""
    output = tmp_path / "out.h5"
    first = run_cmuts(data, output)

    second = run_cmuts(data, output, overwrite=True, min_mapq=60)
    alone = run_cmuts(data, tmp_path / "alone.h5", min_mapq=60)

    assert second != first, "the two runs under test count the same"
    assert outputs_agree(output, tmp_path / "alone.h5")


def test_a_run_that_would_fail_destroys_nothing(data, datasets, tmp_path):
    output = tmp_path / "out.h5"
    run_cmuts(data, output)
    before = output.read_bytes()

    attempt = try_cmuts(datasets("sparse"), output)

    assert attempt.returncode != 0
    assert output.read_bytes() == before


def test_an_empty_file_is_replaced_without_overwrite(data, tmp_path):
    """An empty file is what mktemp and shell redirection leave behind."""
    output = tmp_path / "reserved.h5"
    output.touch()

    assert run_cmuts(data, output).kept > 0


def test_a_file_that_is_not_an_output_is_left_intact(data, tmp_path):
    notes = tmp_path / "notes.txt"
    notes.write_text("months of irreplaceable notes\n")

    attempt = try_cmuts(data, notes)

    assert attempt.returncode != 0
    assert notes.read_text() == "months of irreplaceable notes\n"
