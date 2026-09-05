"""Serving a report over outputs.

Smoke tests only: each starts cmuts plot over some inputs and checks that a
report is served, or that a bad invocation is refused before serving. What the
page holds is not tested here.
"""

from __future__ import annotations

import contextlib
import json
import os
import re
import subprocess
import sys
import urllib.request
from urllib.parse import quote

import pytest

from inputs import not_hdf5, random_fields
from programs import CMUTS_PLOT

# How long to wait for the server to print its address, and for each request.
STARTUP_SECONDS = 30
REQUEST_SECONDS = 30


def script_environment() -> dict:
    """The environment the script runs under, with this interpreter's
    directory first on PATH. The script's /usr/bin/env python3 then resolves
    to the Python holding the test dependencies."""
    environment = os.environ.copy()
    environment["PATH"] = os.path.dirname(sys.executable) + os.pathsep + environment["PATH"]

    return environment


def _address_of(process) -> str:
    """Reads the server's address from what it prints on startup."""
    for _ in range(100):
        line = process.stdout.readline()

        if not line:
            raise AssertionError("the server exited before printing its address")

        found = re.search(r"http://\S+", line)

        if found:
            return found.group(0).rstrip("/")

    raise AssertionError("the server printed no address")


@contextlib.contextmanager
def serving(*arguments):
    """Runs cmuts plot over the arguments, yielding the address it serves at."""
    process = subprocess.Popen(
        [*CMUTS_PLOT, "--port", "0", *map(str, arguments)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        env=script_environment(),
    )

    try:
        yield _address_of(process)
    finally:
        process.terminate()
        process.wait(timeout=STARTUP_SECONDS)


def fetch(address, path=""):
    """One request to the server: the status and the body."""
    with urllib.request.urlopen(address + path, timeout=REQUEST_SECONDS) as reply:
        return reply.status, reply.read()


def assert_serves_a_report(address):
    """Asserts that the address answers with a page and a healthy status."""
    status, body = fetch(address)

    assert status == 200
    assert body

    assert fetch(address, "/health")[0] == 200


def spec_path(spec: dict) -> str:
    return "/figure?spec=" + quote(json.dumps(spec))


def attempt(arguments):
    """Runs an invocation that must exit on its own, bounded so a run that
    wrongly serves fails rather than hangs."""
    return subprocess.run([*CMUTS_PLOT, *map(str, arguments)],
                          capture_output=True, text=True, timeout=STARTUP_SECONDS,
                          env=script_environment())


# ---------------------------------------------------------------------------
# Inputs that are served
# ---------------------------------------------------------------------------


def test_one_file_is_served(build):
    with serving(build()) as address:
        assert_serves_a_report(address)


def test_several_labeled_files_are_served(build):
    files = [build(random_fields(seed=60 + i)) for i in range(3)]
    labels = []

    for name in ("Treated", "Untreated", "Subtracted"):
        labels += ["--label", name]

    with serving(*files, *labels) as address:
        assert_serves_a_report(address)


def test_files_of_different_shapes_are_served(build):
    small = build(n_refs=1, cap=4)
    large = build(n_refs=7, cap=9)

    with serving(small, large) as address:
        assert_serves_a_report(address)


def test_a_figure_is_served(build):
    """One end-to-end figure request, so the figure path is exercised without
    asserting on what it draws."""
    with serving(build()) as address:
        status, body = fetch(address, spec_path(
            {"kind": "profile", "cond": "input0", "ref": 0}))

        assert status == 200
        assert body


# ---------------------------------------------------------------------------
# Invocations that are refused
# ---------------------------------------------------------------------------


def test_a_file_that_is_not_an_output_is_refused(tmp_path):
    failed = attempt([not_hdf5(tmp_path)])

    assert failed.returncode != 0


def test_a_missing_file_is_refused(tmp_path):
    failed = attempt([tmp_path / "absent.h5"])

    assert failed.returncode != 0


def test_more_labels_than_files_is_refused(build):
    failed = attempt([build(), "--label", "a", "--label", "b"])

    assert failed.returncode != 0


def test_a_file_is_required(tmp_path):
    given = attempt([])

    assert given.returncode == 2


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("flag", ["--version", "--help"])
def test_the_informational_flags_exit_cleanly(flag):
    given = attempt([flag])

    assert given.returncode == 0
    assert given.stdout
