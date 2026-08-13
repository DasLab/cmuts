"""What happens to a file already at the output path.

Replacing a previous result requires --overwrite, a run costing far more than
the command that starts it.
"""

import pytest

from datasets import DATASETS
from support import outputs_agree, run_cmuts, try_cmuts


def other_than(name):
    """The next dataset in the catalogue, so that a failing run is given
    something other than what wrote the file it must not touch."""
    names = sorted(DATASETS)

    return names[(names.index(name) + 1) % len(names)]


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_an_existing_output_is_not_replaced_without_overwrite(datasets, tmp_path, name):
    data = datasets(name)
    output = tmp_path / "out.h5"
    run_cmuts(data, output)
    before = output.read_bytes()

    attempt = try_cmuts(data, output)

    assert attempt.returncode != 0
    assert "already holds data" in attempt.stderr
    assert output.read_bytes() == before, "the first result is untouched"


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_overwrite_replaces_an_existing_output(datasets, checked, tmp_path, name):
    """The two runs are given different criteria, so that agreeing with one
    identifies which result was left at the path. Where the criteria make no
    difference to what a dataset counts, agreement identifies nothing."""
    data = datasets(name)
    output = tmp_path / "out.h5"
    first = run_cmuts(data, output)

    second = run_cmuts(data, output, overwrite=True, min_mapq=60)
    run_cmuts(data, tmp_path / "separate.h5", min_mapq=60)

    checked(second != first)

    assert outputs_agree(output, tmp_path / "separate.h5")


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_run_that_would_fail_destroys_nothing(datasets, tmp_path, name):
    """The second run is given another dataset, so a result written before the
    refusal would be the wrong shape as well as the wrong values."""
    output = tmp_path / "out.h5"
    run_cmuts(datasets(name), output)
    before = output.read_bytes()

    attempt = try_cmuts(datasets(other_than(name)), output)

    assert attempt.returncode != 0
    assert output.read_bytes() == before


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_an_empty_file_is_replaced_without_overwrite(datasets, tmp_path, name):
    data = datasets(name)
    output = tmp_path / "reserved.h5"
    output.touch()

    # Not what it counted: a dataset whose reads are all rejected still leaves
    # a row for every reference any of them reached.
    assert run_cmuts(data, output).rows == data.touched


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_file_that_is_not_an_output_is_left_intact(datasets, tmp_path, name):
    data = datasets(name)
    notes = tmp_path / "notes.txt"
    notes.write_text("months of irreplaceable notes\n")

    attempt = try_cmuts(data, notes)

    assert attempt.returncode != 0
    assert notes.read_text() == "months of irreplaceable notes\n"
