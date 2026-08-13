"""Shared fixtures, and the check that every test asserts something somewhere.

Datasets are generated once per session and cached by name, since building one
costs more than every test that reads it.

Every test runs over every dataset. A test states a contract; the datasets are
the range of inputs it is tried on. No dataset holds every case, so a test is
often vacuous on one of them: there is nothing there to assert over. It records
that through the checked fixture instead of failing. A test vacuous on every
dataset asserts nothing at all, and fails the run.
"""

import shutil
from collections import defaultdict

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
        pytest.exit("samtools is required; it is what the tests check against", returncode=2)


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


# How many cases each test asserted over, keyed by test and not by parameter,
# so that every dataset a test ran over is judged together.
_CHECKED = defaultdict(list)


def _amount(cases) -> int:
    """Returns how many cases there are: the size of a container, or the number
    itself.

    A boolean mask is neither, and raises. Its length is the same whether it
    selects everything or nothing, so len() would report full coverage for a
    test that asserted over nothing. Counting its True values instead would be
    wrong the other way: a case need not be a nonzero value, and
    test_weights.py asserts that a set of rates are all zero.
    """
    if getattr(cases, "dtype", None) == bool:
        raise TypeError("checked takes the cases or how many, not a mask over them")

    return len(cases) if hasattr(cases, "__len__") else int(cases)


@pytest.fixture
def checked(request):
    """Records how many cases a test asserts over, and returns them.

    Takes the place of asserting that the case appears in the data:
    `for name in checked(missing)` runs over the references no read reached,
    and zero of them marks the test vacuous on this dataset.
    """
    name = request.node.originalname or request.node.name

    def record(cases):
        _CHECKED[name].append(_amount(cases))
        return cases

    return record


def _vacuous_everywhere():
    """Returns the tests that were vacuous on every dataset.

    Only a run over the whole catalogue shows this: a shorter one may have left
    out the datasets holding the case.
    """
    return sorted(
        name for name, counts in _CHECKED.items()
        if len(counts) >= len(DATASETS) and not any(counts)
    )


def pytest_terminal_summary(terminalreporter):
    vacuous = _vacuous_everywhere()

    if vacuous:
        terminalreporter.section("Vacuous on every dataset")
        for name in vacuous:
            terminalreporter.line(name)


def pytest_sessionfinish(session, exitstatus):
    # A run that already failed says so, and one that stopped early is no
    # evidence either way: its tests did not all run.
    if exitstatus == 0 and _vacuous_everywhere():
        session.exitstatus = 1
