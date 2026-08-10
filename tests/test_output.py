"""What happens to whatever is already at the output path.

A run costs far more than the command that starts it, so replacing a previous
result is something to ask for rather than something to discover afterwards.
"""

import pytest

from support import outputs_agree, run_cmuts, try_cmuts


def test_refuses_to_replace_an_existing_file(data, tmp_path):
    output = tmp_path / "out.h5"
    run_cmuts(data, output)
    before = output.read_bytes()

    attempt = try_cmuts(data, output)

    assert attempt.returncode != 0
    assert "already holds data" in attempt.stderr
    assert output.read_bytes() == before, "the first result is untouched"


def test_overwrite_replaces_it(data, tmp_path):
    """Not that the second run succeeds, but that what is left at the path is
    its result and not the first's. The two are made to differ, or agreeing
    with either would say nothing."""
    output = tmp_path / "out.h5"
    first = run_cmuts(data, output)

    second = run_cmuts(data, output, overwrite=True, min_mapq=60)
    alone = run_cmuts(data, tmp_path / "alone.h5", min_mapq=60)

    assert second != first, "the two runs under test count the same"
    assert outputs_agree(output, tmp_path / "alone.h5")


def test_a_run_that_would_fail_destroys_nothing(data, datasets, tmp_path):
    """The mistake the reference check exists to catch used to truncate the
    output before the check ever ran."""
    output = tmp_path / "out.h5"
    run_cmuts(data, output)
    before = output.read_bytes()

    attempt = try_cmuts(datasets("sparse"), output)

    assert attempt.returncode != 0
    assert output.read_bytes() == before


def test_an_empty_file_is_not_worth_protecting(data, tmp_path):
    """mktemp and shell redirection both leave one behind, and there is nothing
    in it to lose."""
    output = tmp_path / "reserved.h5"
    output.touch()

    assert run_cmuts(data, output).kept > 0


def test_a_file_that_is_not_ours_is_left_alone(data, tmp_path):
    """The path need not hold an output, or anything readable at all."""
    notes = tmp_path / "notes.txt"
    notes.write_text("months of irreplaceable notes\n")

    attempt = try_cmuts(data, notes)

    assert attempt.returncode != 0
    assert notes.read_text() == "months of irreplaceable notes\n"
