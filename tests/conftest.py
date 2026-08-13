"""Shared fixtures and vacuity checks.

Datasets are generated once per session and cached by name. Tests which parametrize over datasets must make use of the falsifiable fixture to ensure at least one dataset would detect an issue in the contract the test checks for. A test vacuous on every dataset fails the run rather than silently passing.
"""

import shutil
from collections import defaultdict

import numpy as np
import pytest

# support carries assertions of its own, and pytest rewrites them only where it
# is told to before the module is first imported.
pytest.register_assert_rewrite("support")

from datasets import DATASETS  # noqa: E402
from support import PROGRAMS, ROOT, generate, located  # noqa: E402


def pytest_configure(config):
    missing = [name for name in PROGRAMS if located(name) is None]
    if missing:
        pytest.exit(
            f"run make check, or put a build directory under {ROOT} first on PATH: "
            f"{', '.join(missing)} not found there",
            returncode=2,
        )

    if not shutil.which("samtools"):
        pytest.exit("samtools is required", returncode=2)


@pytest.fixture(scope="session")
def datasets(tmp_path_factory):
    """Looks a dataset up by name, building it the first time it is asked for."""
    directory = tmp_path_factory.mktemp("datasets")
    built = {}

    def get(name):
        if name not in built:
            built[name] = generate(directory, name, **DATASETS[name])
        return built[name]

    return get


# Whether each test had anything to assert over, keyed by test and not by
# parameter, so that every dataset a test ran over is judged together.
_DECLARED = defaultdict(list)


def _named(node) -> str:
    """Returns the test's name without the parameters it was called with."""
    return node.originalname or node.name


def _dataset_under_test(node):
    """Returns the dataset a test was called with, or None where it ran over
    none. A test runs over a dataset when its name parameter is one in the
    catalogue."""
    callspec = getattr(node, "callspec", None)
    name = callspec.params.get("name") if callspec else None

    return name if name in DATASETS else None


@pytest.fixture(autouse=True)
def _declaration(request):
    """Records a test that ends without declaring as vacuous on this dataset.

    A missing declaration counts as False, so a test that never calls
    falsifiable fails the run in the same way as one that declares False on
    every dataset.
    """
    if _dataset_under_test(request.node) is None:
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
    `falsifiable(len(missing) > 0)` the cases are the references no read
    reached; where a dataset holds none, the test is vacuous on it.
    """
    name = _named(request.node)

    def record(has_cases):
        if not isinstance(has_cases, (bool, np.bool_)):
            raise TypeError("falsifiable takes whether there are cases, not the cases")

        _DECLARED[name].append(bool(has_cases))

    return record


def _vacuous_everywhere():
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
    # Failures take precedence over vacuity
    if exitstatus == 0 and _vacuous_everywhere():
        session.exitstatus = 1
