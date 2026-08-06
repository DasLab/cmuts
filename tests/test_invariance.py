"""Results depend on the input, not on how the work was divided."""

import pytest

from support import outputs_agree, run_cmuts


@pytest.mark.parametrize("shape", ["plain", "sparse", "single"])
def test_worker_count_does_not_change_the_result(datasets, tmp_path, shape):
    data = datasets(shape)

    run_cmuts(data, tmp_path / "one.h5", workers=1, min_mapq=10)
    run_cmuts(data, tmp_path / "many.h5", workers=16, min_mapq=10)

    assert outputs_agree(tmp_path / "one.h5", tmp_path / "many.h5")


@pytest.mark.parametrize("batch", [1, 7, 4096])
def test_batch_size_does_not_change_the_result(datasets, tmp_path, batch):
    """Reads and the carriers holding them move in batches; at a batch of one
    the batched and single-item paths must coincide exactly."""
    data = datasets("plain")

    run_cmuts(data, tmp_path / "default.h5", workers=4)
    run_cmuts(data, tmp_path / "sized.h5", workers=4, batch=batch)

    assert outputs_agree(tmp_path / "default.h5", tmp_path / "sized.h5")


def test_a_queue_of_one_does_not_change_the_result(datasets, tmp_path):
    """Forces the loader to block on the carrier pool for every read."""
    data = datasets("plain")

    run_cmuts(data, tmp_path / "default.h5", workers=4)
    run_cmuts(data, tmp_path / "starved.h5", workers=4, queue_capacity=1, batch=1)

    assert outputs_agree(tmp_path / "default.h5", tmp_path / "starved.h5")
