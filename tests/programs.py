"""Tools for locating and running programs.

Every subprocess the suite starts is started here.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from outputs import read_summary

ROOT = Path(__file__).resolve().parent.parent

CMUTS_HMM = "cmuts-hmm"
CMUTS_GEN = "cmuts-gen"
CMUTS_SUB = "cmuts-sub"
PROGRAMS = (CMUTS_HMM, CMUTS_GEN, CMUTS_SUB)

DEFAULT_WORKERS = 4


def locate(name: str) -> Path | None:
    """Finds a program on PATH, disregarding any copy from outside this
    repository.

    make check puts the build directory first on PATH, which decides
    between an ordinary build and a sanitized one.
    """
    found = shutil.which(name)

    if found is None:
        return None

    path = Path(found).resolve()

    return path if ROOT in path.parents else None


def _words(command) -> list:
    return [str(word) for word in command]


def execute(command):
    """Runs a command, raising where it fails."""
    return subprocess.run(_words(command), check=True, capture_output=True, text=True)


def attempt(command):
    """Runs a command whether or not it succeeds."""
    return subprocess.run(_words(command), capture_output=True, text=True)


def execute_into(path, command):
    """Runs a command with its standard output written to a file."""
    with open(path, "wb") as handle:
        subprocess.run(_words(command), check=True, stdout=handle,
                       stderr=subprocess.DEVNULL)

    return path


def reported_version(program: str) -> str:
    """Returns the version a program prints for --version, which it gives after
    its own name."""
    return execute([program, "--version"]).stdout.split()[-1]


def samtools(*arguments) -> str:
    """Runs samtools and returns what it printed."""
    return execute(["samtools", *arguments]).stdout


def samtools_into(path, *arguments):
    """Runs samtools with its output written to a file."""
    return execute_into(path, ["samtools", *arguments])


def option(name: str) -> str:
    """Returns the command-line option a keyword argument names."""
    return "--" + name.replace("_", "-")


def _options(given: dict) -> list:
    """Returns the words a set of keyword arguments contributes. A flag has
    no value of its own."""
    words = []

    for name, value in given.items():
        words += [option(name)] if value is True else [option(name), value]

    return words


def _cmuts_command(data, output, workers: int, options: dict) -> list:
    return [CMUTS_HMM, "-f", data.fasta, "-o", output, "-j", workers,
            *_options(options), *data.bams]


def run_cmuts(data, output, workers: int = DEFAULT_WORKERS, **options):
    """Counts an alignment, returning the read counts written to the output."""
    execute(_cmuts_command(data, output, workers, options))

    return read_summary(output)


def try_cmuts(data, output, workers: int = DEFAULT_WORKERS, **options):
    """Counts an alignment whether or not the run succeeds."""
    return attempt(_cmuts_command(data, output, workers, options))


def _subtract_command(treated, untreated, output, options: dict) -> list:
    return [CMUTS_SUB, "-o", output, *_options(options), treated, untreated]


def run_subtract(treated, untreated, output, **options):
    """Subtracts an untreated background from a treated run, returning the path
    the result was written to."""
    execute(_subtract_command(treated, untreated, output, options))

    return output


def try_subtract(treated, untreated, output, **options):
    """Subtracts a background whether or not the run succeeds."""
    return attempt(_subtract_command(treated, untreated, output, options))


def run_generator(prefix, parameters: dict):
    """Writes an alignment and its reference, both named after the prefix."""
    execute([CMUTS_GEN, "-o", prefix, *_options(parameters)])
