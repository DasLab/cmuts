"""The pair HMM's rates, read from a file.

`--params` replaces the built-in rates. A file is refused where a rate is out of
range or the two opening rates sum above one. Anything else is read, and whether
a read can be scored against the model it describes is settled per read: a model
that forbids what a read contains gives that read a probability of zero, and the
read is counted as rejected.
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

# A model that forbids indels, by never opening one. A read whose indels are short
# enough for the band to cross is scored as a run of mismatches; a read carrying a
# longer one has no path through the model.
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
    them back is a run given no overrides."""
    dumped = execute_into(tmp_path / "dumped.txt", [*CMUTS_HMM, "--dump-params"])

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
    """Refused before the alignments are opened, so no output file is written."""
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
def test_a_read_with_no_path_is_not_counted(data, falsifiable, tmp_path):
    params = write_params(tmp_path / "params.txt", **NO_SECOND_BASE)

    summary = read_summary(run_cmuts(data, tmp_path / "out.h5", params=params))

    falsifiable(data.mapped > 0)

    assert summary.kept == 0


@pytest.mark.parametrize("fmt", [NATIVE], indirect=True)
def test_a_read_with_no_path_is_counted_as_rejected(data, falsifiable, tmp_path):
    """Every mapped read lands in one of the two counts whether or not the model
    can score it, so a run against rates that fit no read still accounts for the
    file it read."""
    params = write_params(tmp_path / "params.txt", **NO_SECOND_BASE)

    summary = read_summary(run_cmuts(data, tmp_path / "out.h5", params=params))

    falsifiable(data.mapped > 0)

    assert summary.kept + summary.rejected == data.mapped


# The bases each deletion spans. The band is the other half of what makes a path,
# so a band reaching further than this can cross one and a narrower band cannot.
DELETION_LENGTH = 3


@pytest.fixture
def short_deletions(tmp_path):
    """Returns an alignment whose reads carry short deletions and no insertions."""
    return generate(tmp_path, "short_deletions", seed=7, references=3, ref_length=300,
                    reads_per_ref=40, insertions=0, deletions=1,
                    deletion_length=DELETION_LENGTH)


def scored(data, path, **options) -> int:
    """The reads a run counted."""
    return read_summary(run_cmuts(data, path, **options)).kept


def test_a_narrow_band_scores_fewer_reads(short_deletions, tmp_path):
    """A band that cannot reach across a deletion leaves the reads carrying one no
    path, and a band that can scores them as mismatches instead."""
    params = write_params(tmp_path / "params.txt", **NO_INDELS)

    narrow = scored(short_deletions, tmp_path / "narrow.h5", params=params,
                    band=str(DELETION_LENGTH - 1))
    wide = scored(short_deletions, tmp_path / "wide.h5", params=params,
                  band=str(DELETION_LENGTH + 1))

    assert narrow < wide


def test_forbidding_indels_scores_fewer_reads(catalogue, tmp_path):
    """Compared against the rates a default run uses, over the one dataset, so that
    what the two disagree on is the model and not the reads."""
    data = catalogue("indels", NATIVE)
    params = write_params(tmp_path / "params.txt", **NO_INDELS)

    forbidden = scored(data, tmp_path / "forbidden.h5", params=params)
    allowed = scored(data, tmp_path / "allowed.h5")

    assert forbidden < allowed


def test_forbidding_indels_scores_a_run_without_them(catalogue, tmp_path):
    """The same rates that leave most of a dataset carrying indels unscored: whether
    a model can score a run's reads is a property of the two together."""
    params = write_params(tmp_path / "params.txt", **NO_INDELS)

    assert scored(catalogue("clean", NATIVE), tmp_path / "out.h5", params=params) > 0
