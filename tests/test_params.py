"""The pair HMM's rates, read from a file.

`--params` replaces the built-in rates. A file is refused where a rate is out of
range or the two opening rates sum above one. Anything else is read, and whether
a read can be scored against the model it describes is settled per read: a model
that forbids what a read contains gives that read a probability of zero, which
stops the run.
"""

from __future__ import annotations

import pytest

from alignments import NATIVE, generate
from outputs import outputs_agree, read_summary
from programs import CMUTS_HMM, execute_into, run_cmuts, try_cmuts

# Every key the reader knows, at the rates a default run uses.
USABLE = {
    "deletion-open":     "0.001",
    "deletion-extend":   "0.35",
    "insertion-open":    "0.0002",
    "insertion-extend":  "0.35",
    "modification-rate": "0.005",
}

# A model no read of more than one base can be scored against: every step out of
# a match is a deletion, which consumes a reference base and not a read one, so
# there is no way to reach the second base of a read. Legal to the reader, every
# rate being in range and the openings summing to one.
NO_SECOND_BASE = {"deletion-open": "1.0", "insertion-open": "0.0"}

# A model that forbids indels, by never opening one. A read carrying an indel has
# no path through it; a read without one is scored as usual.
NO_INDELS = {"deletion-open": "0", "insertion-open": "0"}

# Files the reader itself refuses.
REFUSED = {
    "above one":    {"modification-rate": "1.5"},
    "below zero":   {"deletion-extend": "-0.1"},
    "openings sum": {"deletion-open": "0.7", "insertion-open": "0.7"},
    "unknown key":  {"substitution-rate": "0.1"},
    "not a number": {"deletion-open": "often"},
}


def write_params(path, **changed):
    """Writes a parameter file, starting from the rates a default run uses."""
    rates = USABLE | changed
    path.write_text("".join(f"{name} {rate}\n" for name, rate in rates.items()))

    return path


# ---------------------------------------------------------------------------
# The rates a default run uses
# ---------------------------------------------------------------------------


def test_dumped_defaults_reproduce_a_default_run(data, falsifiable, tmp_path):
    """--dump-params writes the rates in the form --params reads, so a run given
    them back is a run given nothing."""
    dumped = execute_into(tmp_path / "dumped.txt", [CMUTS_HMM, "--dump-params"])

    run_cmuts(data, tmp_path / "plain.h5")
    run_cmuts(data, tmp_path / "given.h5", params=dumped)

    falsifiable(data.mapped > 0)

    assert outputs_agree(tmp_path / "plain.h5", tmp_path / "given.h5")


def test_the_rates_are_read_in_any_order(data, falsifiable, tmp_path):
    forward = write_params(tmp_path / "forward.txt")
    backward = tmp_path / "backward.txt"
    lines = forward.read_text().splitlines(keepends=True)
    backward.write_text("".join(reversed(lines)))

    run_cmuts(data, tmp_path / "forward.h5", params=forward)
    run_cmuts(data, tmp_path / "backward.h5", params=backward)

    falsifiable(data.mapped > 0)

    assert outputs_agree(tmp_path / "forward.h5", tmp_path / "backward.h5")


# ---------------------------------------------------------------------------
# Files the reader refuses
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fmt", [NATIVE], indirect=True)
@pytest.mark.parametrize("wrong", sorted(REFUSED))
def test_a_refused_file_leaves_no_output(data, falsifiable, tmp_path, wrong):
    """Refused before the alignments are opened, so nothing is written whatever
    the dataset holds."""
    params = write_params(tmp_path / "params.txt", **REFUSED[wrong])
    output = tmp_path / "out.h5"

    failed = try_cmuts(data, output, params=params)

    falsifiable(True)

    assert failed.returncode != 0
    assert not output.exists()


def test_a_missing_file_is_refused(data, falsifiable, tmp_path):
    failed = try_cmuts(data, tmp_path / "out.h5", params=tmp_path / "absent.txt")

    falsifiable(True)

    assert failed.returncode != 0


# ---------------------------------------------------------------------------
# Reads the model gives no path
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fmt", [NATIVE], indirect=True)
def test_a_read_with_no_path_stops_the_run(data, falsifiable, tmp_path):
    """A dataset from which no read is counted has nothing to score and runs to
    the end, which is what makes this conditional rather than unconditional."""
    params = write_params(tmp_path / "params.txt", **NO_SECOND_BASE)
    counted = read_summary(run_cmuts(data, tmp_path / "plain.h5")).kept

    failed = try_cmuts(data, tmp_path / "out.h5", params=params)

    falsifiable(counted > 0)

    if counted > 0:
        assert failed.returncode != 0
    else:
        assert failed.returncode == 0


# The bases each deletion spans. The band is the other half of what makes a path,
# so a band reaching further than this can cross one and a narrower band cannot.
DELETION_LENGTH = 3


@pytest.fixture
def short_deletions(tmp_path):
    """Returns an alignment whose reads carry short deletions and no insertions."""
    return generate(tmp_path, "short_deletions", seed=7, references=3, ref_length=300,
                    reads_per_ref=40, insertions=0, deletions=1,
                    deletion_length=DELETION_LENGTH)


def test_a_narrow_band_gives_no_path(short_deletions, tmp_path):
    params = write_params(tmp_path / "params.txt", **NO_INDELS)

    failed = try_cmuts(short_deletions, tmp_path / "out.h5", params=params,
                       band=str(DELETION_LENGTH - 1))

    assert failed.returncode != 0


def test_a_wide_band_admits_a_path(short_deletions, tmp_path):
    """The deletion is scored as mismatches instead, so the rates that stop the
    narrower band run to the end."""
    params = write_params(tmp_path / "params.txt", **NO_INDELS)

    run = try_cmuts(short_deletions, tmp_path / "out.h5", params=params,
                    band=str(DELETION_LENGTH + 1))

    assert run.returncode == 0


def test_forbidding_indels_stops_a_run_with_them(catalogue, tmp_path):
    params = write_params(tmp_path / "params.txt", **NO_INDELS)

    failed = try_cmuts(catalogue("indels", NATIVE), tmp_path / "out.h5", params=params)

    assert failed.returncode != 0


def test_forbidding_indels_allows_a_run_without_them(catalogue, tmp_path):
    """The same rates that stop a dataset carrying indels: whether a model can
    score a run's reads is a property of the two together."""
    params = write_params(tmp_path / "params.txt", **NO_INDELS)

    run = try_cmuts(catalogue("clean", NATIVE), tmp_path / "out.h5", params=params)

    assert run.returncode == 0
