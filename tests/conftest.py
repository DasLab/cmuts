"""Shared fixtures.

Datasets are generated once per session and cached by name, since building one
costs more than every test that reads it.
"""

import shutil

import pytest

# support carries assertions of its own, and pytest rewrites them only where it
# is told to before the module is first imported.
pytest.register_assert_rewrite("support")

from support import CMUTS, CMUTS_GEN, CMUTS_SUB, generate  # noqa: E402

# Each dataset is one that is easy to get wrong: many references barely covered,
# one reference deeply covered, lengths ranging sixtyfold, reads whose stored
# length far exceeds the span they align to through soft clipping and again
# through long insertions, a file where every read falls below any useful
# threshold, one with no differences from the reference at all, one whose reads
# run past twice the length of the reference they are placed on, and one that is
# both ragged and sparsely covered, so that a reference with padding and no
# reads at all is among them.
DATASETS = {
    "plain":   dict(seed=101, references=40, ref_length="150:600", reads_per_ref=25),
    "sparse":  dict(seed=102, references=800, covered=0.3, reads_per_ref="1:6",
                    ref_length="200:300"),
    "single":  dict(seed=103, references=1, ref_length=600, reads_per_ref=3000),
    "ragged":  dict(seed=104, references=30, ref_length="60:4000", reads_per_ref="10:40"),
    "clipped": dict(seed=105, references=25, reads_per_ref=20, soft_clips=2,
                    soft_clip_length="20:80"),
    "indels":  dict(seed=106, references=25, reads_per_ref=20, insertions="1:3",
                    insertion_length="20:200", deletions="1:2"),
    "lowqual": dict(seed=107, references=15, reads_per_ref=20, mapq=0),
    "clean":   dict(seed=108, references=15, reads_per_ref=20, mismatch_rate=0,
                    insertions=0, deletions=0, soft_clips=0),
    "overflowing": dict(seed=109, references=12, ref_length=60, reads_per_ref=20,
                        read_length="40:60", soft_clips="0:2",
                        soft_clip_length="20:150"),
    "patchy":  dict(seed=110, references=60, ref_length="60:900", covered=0.4,
                    reads_per_ref="5:20"),
    "flat":    dict(seed=111, references=120, ref_length=300, covered=0.4,
                    reads_per_ref="5:20"),
}


def pytest_configure(config):
    missing = [path.name for path in (CMUTS, CMUTS_GEN, CMUTS_SUB) if not path.exists()]
    if missing:
        pytest.exit(f"run make first: {', '.join(missing)} not built", returncode=2)

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
