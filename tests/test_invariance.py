"""The result depends on the input and not on how the work was divided.

Each test compares two runs for equality, so a dataset that produces no rate at
all is still a valid case. A run that kept no reads is not, both runs then
writing empty arrays, so each test declares on what its own run kept.
"""

import pytest

from outputs import outputs_agree
from programs import run_cmuts

# The divisions of work a run is tried under.
DECODE_THREADS = [0, 1, 8]
BATCH_SIZES = [1, 7, 4096]


@pytest.fixture(scope="module")
def baseline(tmp_path_factory):
    """Returns a function giving the result each variation is compared against,
    run once per dataset."""
    directory = tmp_path_factory.mktemp("baseline")
    made = {}

    def get(data):
        if data.name not in made:
            output = directory / f"{data.name}.h5"
            run_cmuts(data, output)
            made[data.name] = output

        return made[data.name]

    return get


def test_worker_count_does_not_change_the_result(data, falsifiable, tmp_path):
    one = run_cmuts(data, tmp_path / "one.h5", workers=1, min_mapq=10)
    run_cmuts(data, tmp_path / "many.h5", workers=16, min_mapq=10)

    falsifiable(one.kept > 0)

    assert outputs_agree(tmp_path / "one.h5", tmp_path / "many.h5")


@pytest.mark.parametrize("threads", DECODE_THREADS)
def test_decode_threads_do_not_change_the_result(data, falsifiable, baseline, tmp_path,
                                                 threads):
    varied = run_cmuts(data, tmp_path / "threaded.h5", decode_threads=threads)

    falsifiable(varied.kept > 0)

    assert outputs_agree(baseline(data), tmp_path / "threaded.h5")


@pytest.mark.parametrize("batch", BATCH_SIZES)
def test_batch_size_does_not_change_the_result(data, falsifiable, baseline, tmp_path,
                                               batch):
    varied = run_cmuts(data, tmp_path / "sized.h5", batch=batch)

    falsifiable(varied.kept > 0)

    assert outputs_agree(baseline(data), tmp_path / "sized.h5")


def test_a_queue_of_one_does_not_change_the_result(data, falsifiable, baseline,
                                                   tmp_path):
    varied = run_cmuts(data, tmp_path / "starved.h5", queue_capacity=1, batch=1)

    falsifiable(varied.kept > 0)

    assert outputs_agree(baseline(data), tmp_path / "starved.h5")


def test_one_reference_in_flight_does_not_change_the_result(data, falsifiable, baseline,
                                                            tmp_path):
    varied = run_cmuts(data, tmp_path / "single-file.h5", live_refs=1)

    falsifiable(varied.kept > 0)

    assert outputs_agree(baseline(data), tmp_path / "single-file.h5")
