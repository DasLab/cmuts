"""Shared fixtures and the check for vacuous tests.

Datasets are generated once per session and cached by name. A test that runs
over the catalogue must call the falsifiable fixture to state whether the
dataset it was given can detect a violation of the contract under test. A test
that is vacuous on every dataset fails the run.
"""

import itertools
import shutil
from collections import defaultdict

import numpy as np
import pytest

# oracle contains assertions of its own, and pytest rewrites them only where it
# is told to before the module is first imported.
pytest.register_assert_rewrite("oracle")

from alignments import FORMATS, NATIVE, convert_format, generate  # noqa: E402
from datasets import DATASETS  # noqa: E402
from inputs import CAP, N_REFS  # noqa: E402
from outputs import write_output  # noqa: E402
from programs import PROGRAMS, ROOT, locate  # noqa: E402


def pytest_configure(config):
    missing = [name for name in PROGRAMS if locate(name) is None]
    if missing:
        pytest.exit(
            f"run make check, or put a build directory under {ROOT} first on PATH: "
            f"{', '.join(missing)} not found there",
            returncode=2,
        )

    if not shutil.which("samtools"):
        pytest.exit("samtools is required", returncode=2)


@pytest.fixture(scope="session")
def catalogue(tmp_path_factory):
    """Returns a function that looks a dataset up by name and format,
    generating and converting it on first use."""
    directory = tmp_path_factory.mktemp("datasets")
    built = {}

    def get(name, fmt):
        if (name, fmt) not in built:
            native = built.setdefault(
                (name, NATIVE), generate(directory, name, **DATASETS[name]))
            built[(name, fmt)] = convert_format(native, directory, fmt)

        return built[(name, fmt)]

    return get


@pytest.fixture(params=FORMATS)
def fmt(request):
    """One format htslib reads. A contract holds for every format the alignments
    arrive in, so every test over the catalogue is run against each.

    A test overrides this fixture when it requires specific formats.
    """
    return request.param


@pytest.fixture(params=sorted(DATASETS))
def data(request, catalogue, fmt):
    """One dataset of the catalogue, in one format. A test taking this fixture
    runs over every combination of the two."""
    return catalogue(request.param, fmt)


@pytest.fixture(params=["plain", "chunked"])
def storage(request):
    """The two ways an input may be stored. cmuts hmm writes chunked, shuffled
    and deflated, and the result must be the same either way, so every test
    that reads values runs against both."""
    return request.param


@pytest.fixture
def build(tmp_path, storage):
    """Returns a function that writes an input file. Each file is named
    separately, so one test may build several."""
    written = itertools.count()

    def make(values=None, *, n_refs=N_REFS, cap=CAP, unmapped=0):
        return write_output(
            tmp_path / f"input{next(written)}.h5",
            n_refs=n_refs, cap=cap, values=values, unmapped=unmapped,
            storage=storage,
        )

    return make


# Whether each test had anything to assert over, keyed by test and not by
# parameter, so that every dataset a test ran over is judged together.
_DECLARED = defaultdict(list)


def _named(node) -> str:
    """Returns the test's name without the parameters it was called with."""
    return node.originalname or node.name


def _runs_over_a_dataset(node) -> bool:
    """Returns whether the test was called with a dataset of the catalogue."""
    callspec = getattr(node, "callspec", None)

    return bool(callspec) and callspec.params.get("data") in DATASETS


@pytest.fixture(autouse=True)
def _declaration(request):
    """Records a test that ends without declaring as vacuous on this dataset.

    A missing declaration counts as False, so a test that never calls
    falsifiable fails the run in the same way as one that declares False on
    every dataset.
    """
    if not _runs_over_a_dataset(request.node):
        yield
        return

    name = _named(request.node)
    declared = len(_DECLARED[name])

    yield

    if len(_DECLARED[name]) == declared:
        _DECLARED[name].append(False)


@pytest.fixture
def falsifiable(request):
    """Records whether the test has anything to assert over on this dataset.

    Takes the place of asserting that the case appears in the data. In
    `falsifiable(len(missing) > 0)` the cases are the references that no read
    aligned to; a dataset holding none leaves the test vacuous.
    """
    name = _named(request.node)

    def record(has_cases):
        if not isinstance(has_cases, (bool, np.bool_)):
            raise TypeError("falsifiable takes whether there are cases, not the cases")

        _DECLARED[name].append(bool(has_cases))

    return record


def _vacuous_everywhere() -> list:
    """Returns the tests that were vacuous on every dataset."""
    return sorted(
        name for name, declared in _DECLARED.items()
        if not any(declared)
    )


def pytest_terminal_summary(terminalreporter):
    vacuous = _vacuous_everywhere()

    if vacuous:
        terminalreporter.section("Vacuous on every dataset")
        for name in vacuous:
            terminalreporter.line(name)


def pytest_sessionfinish(session, exitstatus):
    # A failure is reported ahead of any vacuity.
    if exitstatus == 0 and _vacuous_everywhere():
        session.exitstatus = 1
