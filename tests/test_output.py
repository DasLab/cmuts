"""What a run leaves at the output path.

Replacing a previous result requires --overwrite, a run costing far more than
the command that starts it. The program and version stamped on the result are
covered here as well.
"""

from datasets import DATASETS
from outputs import attributes_of, outputs_agree
from programs import CMUTS_HMM, reported_version, run_cmuts, try_cmuts

# Higher than the first run's threshold of zero, so the two runs of an
# overwrite count different reads.
STRICTER_MAPQ = 60

NOTES = "months of irreplaceable notes\n"


def next_dataset(catalogue, data):
    """Returns the dataset after this one in the catalogue, in the same format,
    so that a failing run is given something other than what wrote the file it
    must not touch."""
    names = sorted(DATASETS)
    following = names[(names.index(data.name) + 1) % len(names)]

    return catalogue(following, data.fmt)


def test_an_existing_output_is_not_replaced_without_overwrite(data, falsifiable,
                                                              tmp_path):
    output = tmp_path / "out.h5"
    first = run_cmuts(data, output)
    before = output.read_bytes()

    # The refusal protects a result, so the first run must have written rows.
    falsifiable(first.rows > 0)

    attempt = try_cmuts(data, output)

    assert attempt.returncode != 0
    assert output.read_bytes() == before, "the first result is untouched"


def test_overwrite_replaces_an_existing_output(data, falsifiable, tmp_path):
    """The two runs are given different criteria, so that agreeing with one of
    them identifies which result was left at the path. Where the criteria make
    no difference to what a dataset counts, agreement identifies nothing."""
    output = tmp_path / "out.h5"
    first = run_cmuts(data, output)

    second = run_cmuts(data, output, overwrite=True, min_mapq=STRICTER_MAPQ)
    run_cmuts(data, tmp_path / "separate.h5", min_mapq=STRICTER_MAPQ)

    falsifiable(second != first)

    assert outputs_agree(output, tmp_path / "separate.h5")


def test_a_run_that_would_fail_destroys_nothing(data, falsifiable, catalogue, tmp_path):
    """The second run is given another dataset, so a result written before the
    refusal would be the wrong shape as well as the wrong values."""
    output = tmp_path / "out.h5"
    first = run_cmuts(data, output)
    before = output.read_bytes()

    falsifiable(first.rows > 0)

    attempt = try_cmuts(next_dataset(catalogue, data), output)

    assert attempt.returncode != 0
    assert output.read_bytes() == before


def test_an_empty_file_is_replaced_without_overwrite(data, falsifiable, tmp_path):
    output = tmp_path / "reserved.h5"
    output.touch()

    falsifiable(data.touched > 0)

    # Compared against the references any read reached: a dataset whose reads
    # are all rejected still gets a row for every reference they aligned to.
    assert run_cmuts(data, output).rows == data.touched


def test_a_file_that_is_not_an_output_is_left_intact(data, falsifiable, tmp_path):
    notes = tmp_path / "notes.txt"
    notes.write_text(NOTES)

    # cmuts-hmm refuses on the file already at the path and never reads the
    # alignments, so no dataset can leave this test with nothing to assert.
    falsifiable(True)

    attempt = try_cmuts(data, notes)

    assert attempt.returncode != 0
    assert notes.read_text() == NOTES


def test_outputs_are_labelled_with_the_program(data, falsifiable, tmp_path):
    output = tmp_path / "out.h5"

    # The attribute is written by every run, whatever it counted.
    falsifiable(True)

    run_cmuts(data, output)

    assert attributes_of(output)["program"] == CMUTS_HMM


def test_outputs_are_versioned(data, falsifiable, tmp_path):
    """Compared against the version the program prints, so that a release
    cannot change one without the other."""
    output = tmp_path / "out.h5"

    falsifiable(True)

    run_cmuts(data, output)

    assert attributes_of(output)["version"] == reported_version(CMUTS_HMM)
