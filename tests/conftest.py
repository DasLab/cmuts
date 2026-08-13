"""Shared fixtures.

Datasets are generated once per session and cached by name, since building one
costs more than every test that reads it.
"""

import shutil

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


@pytest.fixture
def data(datasets):
    """The everyday dataset, for a test that needs one but not a particular
    one."""
    return datasets("plain")
