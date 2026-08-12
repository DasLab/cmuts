"""Results depend on the input, not on how the work was divided."""

import pytest

from support import outputs_agree, run_cmuts

# Each stresses a different part of the division of work: plain is the ordinary
# case, sparse spends its time crossing between references, single puts every
# read on one and so puts every worker on one accumulator, and patchy is the
# only one carrying references no read reaches and shorter than the widest,
# which the loader opens and closes as it passes them.
SHAPES = ["plain", "sparse", "single", "patchy"]


@pytest.fixture(scope="module")
def unvaried(datasets, tmp_path_factory):
    """The result each variation is compared against, run once per shape."""
    directory = tmp_path_factory.mktemp("unvaried")
    made = {}

    def get(shape):
        if shape not in made:
            output = directory / f"{shape}.h5"
            run_cmuts(datasets(shape), output, workers=4)
            made[shape] = output

        return made[shape]

    return get


@pytest.mark.parametrize("shape", SHAPES)
def test_worker_count_does_not_change_the_result(datasets, tmp_path, shape):
    data = datasets(shape)

    run_cmuts(data, tmp_path / "one.h5", workers=1, min_mapq=10)
    run_cmuts(data, tmp_path / "many.h5", workers=16, min_mapq=10)

    assert outputs_agree(tmp_path / "one.h5", tmp_path / "many.h5")


@pytest.mark.parametrize("shape", SHAPES)
@pytest.mark.parametrize("threads", [0, 1, 8])
def test_decode_threads_do_not_change_the_result(datasets, unvaried, tmp_path,
                                                 shape, threads):
    run_cmuts(datasets(shape), tmp_path / "threaded.h5",
              workers=4, decode_threads=threads)

    assert outputs_agree(unvaried(shape), tmp_path / "threaded.h5")


@pytest.mark.parametrize("shape", SHAPES)
@pytest.mark.parametrize("batch", [1, 7, 4096])
def test_batch_size_does_not_change_the_result(datasets, unvaried, tmp_path,
                                               shape, batch):
    run_cmuts(datasets(shape), tmp_path / "sized.h5", workers=4, batch=batch)

    assert outputs_agree(unvaried(shape), tmp_path / "sized.h5")


@pytest.mark.parametrize("shape", SHAPES)
def test_a_queue_of_one_does_not_change_the_result(datasets, unvaried, tmp_path,
                                                   shape):
    run_cmuts(datasets(shape), tmp_path / "starved.h5",
              workers=4, queue_capacity=1, batch=1)

    assert outputs_agree(unvaried(shape), tmp_path / "starved.h5")


@pytest.mark.parametrize("shape", SHAPES)
def test_one_reference_in_flight_does_not_change_the_result(datasets, unvaried,
                                                            tmp_path, shape):
    run_cmuts(datasets(shape), tmp_path / "single-file.h5", workers=4, live_refs=1)

    assert outputs_agree(unvaried(shape), tmp_path / "single-file.h5")
