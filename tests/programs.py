"""Tools for locating and running programs.

Every subprocess the suite starts is started here.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from outputs import read_summary

ROOT = Path(__file__).resolve().parent.parent

CMUTS = "cmuts"

CMUTS_HMM = (CMUTS, "hmm")
CMUTS_GEN = (CMUTS, "gen")
CMUTS_SUB = (CMUTS, "sub")
CMUTS_DIV = (CMUTS, "div")
CMUTS_NORM = (CMUTS, "norm")
PROGRAMS = (CMUTS,)


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


class ProgramFailed(subprocess.CalledProcessError):
    """A run that exited non-zero. The base class names the status or the signal
    and holds what the program wrote without showing it, which the message here
    adds."""

    def __str__(self) -> str:
        return f"{super().__str__()}\n{self.stderr}"


def _raise_on_failure(result) -> None:
    """Raises where a finished run exited non-zero."""
    if result.returncode != 0:
        raise ProgramFailed(result.returncode, result.args, result.stdout,
                            result.stderr)


def attempt(command):
    """Runs a command whether or not it succeeds."""
    return subprocess.run(_words(command), capture_output=True, text=True)


def execute(command):
    """Runs a command, raising where it fails."""
    result = attempt(command)

    _raise_on_failure(result)

    return result


def execute_into(path, command):
    """Runs a command with its standard output written to a file. Standard error
    is kept apart from it, so that a failure is still reported with what the
    program wrote."""
    with open(path, "wb") as handle:
        result = subprocess.run(_words(command), stdout=handle,
                                stderr=subprocess.PIPE, text=True)

    _raise_on_failure(result)

    return path


def reported_version(program) -> str:
    """Returns the version a program prints for --version, which it gives after
    its own name."""
    return execute([*program, "--version"]).stdout.split()[-1]


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


def _cmuts_command(data, output, options: dict) -> list:
    return [*CMUTS_HMM, "-f", data.fasta, "-o", output, *_options(options), *data.bams]


def run_cmuts(data, output, **options):
    """Counts an alignment, returning the path the result was written to.

    An option the caller does not give is left to the program, so a run here is the run
    a user gets. Read the counts back with read_summary.
    """
    execute(_cmuts_command(data, output, options))

    return output


def try_cmuts(data, output, **options):
    """Counts an alignment whether or not the run succeeds."""
    return attempt(_cmuts_command(data, output, options))


def _subtract_command(treated, untreated, output, options: dict) -> list:
    return [*CMUTS_SUB, "-o", output, *_options(options), treated, untreated]


def run_subtract(treated, untreated, output, **options):
    """Subtracts an untreated background from a treated run, returning the path
    the result was written to."""
    execute(_subtract_command(treated, untreated, output, options))

    return output


def try_subtract(treated, untreated, output, **options):
    """Subtracts a background whether or not the run succeeds."""
    return attempt(_subtract_command(treated, untreated, output, options))


def _divide_command(rates, control, output, options: dict) -> list:
    return [*CMUTS_DIV, "-o", output, *_options(options), rates, control]


def run_divide(rates, control, output, **options):
    """Divides a run by a denatured control, returning the path the result was
    written to."""
    execute(_divide_command(rates, control, output, options))

    return output


def try_divide(rates, control, output, **options):
    """Divides by a control whether or not the run succeeds."""
    return attempt(_divide_command(rates, control, output, options))


def _normalize_command(inputs, outputs, options: dict) -> list:
    """Builds a run over any number of inputs, each paired with its own output
    by a repeat of --output."""
    paired = []

    for output in outputs:
        paired += ["-o", output]

    return [*CMUTS_NORM, *paired, *_options(options), *inputs]


def run_normalize(inputs, outputs, **options):
    """Normalizes every input against one scale, returning the paths written."""
    execute(_normalize_command(inputs, outputs, options))

    return outputs


def try_normalize(inputs, outputs, **options):
    """Normalizes whether or not the run succeeds."""
    return attempt(_normalize_command(inputs, outputs, options))


def run_generator(prefix, parameters: dict):
    """Writes an alignment and its reference, both named after the prefix."""
    execute([*CMUTS_GEN, "-o", prefix, *_options(parameters)])
