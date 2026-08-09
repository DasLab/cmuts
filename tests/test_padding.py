"""What a value of NaN means in the output.

The per-base arrays are rectangular and the references are not, so every row
runs to the length of the longest. NaN is what tells the difference: a position
past the end of its own reference, and a reference no alignment ever named.

A zero is a measurement and says something else entirely -- that the position
was reached and nothing was counted there. The two must not be confused, since
a rate taken over a zero denominator is a division nothing defines, while a
rate taken over a NaN is a NaN and says so. Anything summing these arrays has
to know which it is looking at, so what each one means is fixed here rather
than left to whatever the writer happened to do.
"""

import h5py
import numpy as np
import pytest

from support import references_with_reads, run_cmuts, sequences

# Higher than any mapping quality a read can carry, so every read is turned
# away and the references themselves are all that is left.
REJECTS_EVERYTHING = 61


def per_base(output):
    """The arrays with a value for every reference position."""
    return {name: output[name][:] for name in output if output[name].ndim == 2}


def per_reference(output):
    """The arrays with one value for each reference."""
    return {
        name: output[name][:]
        for name in output
        if output[name].ndim == 1 and np.issubdtype(output[name].dtype, np.floating)
    }


def rows_by_name(output):
    """Which row of the output each reference was written to."""
    names = [name.decode() if isinstance(name, bytes) else str(name)
             for name in output["reference"][:]]
    return {name: i for i, name in enumerate(names)}


@pytest.fixture
def ragged(datasets, tmp_path):
    """Lengths ranging sixtyfold, so most rows are mostly padding."""
    data = datasets("ragged")
    output = tmp_path / "ragged.h5"
    run_cmuts(data, output)
    return data, output


# ---------------------------------------------------------------------------
# A reference shorter than the row it is written to
# ---------------------------------------------------------------------------


def test_positions_past_a_reference_are_nan(ragged):
    data, output = ragged
    lengths = {name: len(seq) for name, seq in sequences(data.fasta).items()}

    with h5py.File(output, "r") as handle:
        row_of = rows_by_name(handle)
        arrays = per_base(handle)
        width = next(iter(arrays.values())).shape[1]

        shorter = [name for name, n in lengths.items() if n < width]
        assert shorter, "the shape under test has no reference shorter than the widest"

        for field, values in arrays.items():
            for name in shorter:
                tail = values[row_of[name]][lengths[name]:]
                assert np.isnan(tail).all(), f"{field}: {name} is padded with {tail[:4]}"


def test_positions_within_a_reference_are_never_nan(ragged):
    """The complement of the above: padding is the only thing NaN marks on a
    reference some read reached, so a NaN inside one is a hole in the counts."""
    data, output = ragged
    lengths = {name: len(seq) for name, seq in sequences(data.fasta).items()}
    reached = references_with_reads(data.bam)

    with h5py.File(output, "r") as handle:
        row_of = rows_by_name(handle)

        for field, values in per_base(handle).items():
            for name in reached:
                within = values[row_of[name]][:lengths[name]]
                assert not np.isnan(within).any(), f"{field}: {name} has a hole in it"


# ---------------------------------------------------------------------------
# A reference no alignment named
# ---------------------------------------------------------------------------


def test_a_reference_no_read_named_is_nan_throughout(datasets, tmp_path):
    """Thirty per cent covered, so most references are named by nothing."""
    data = datasets("sparse")
    output = tmp_path / "sparse.h5"
    run_cmuts(data, output)

    reached = references_with_reads(data.bam)

    with h5py.File(output, "r") as handle:
        row_of = rows_by_name(handle)
        missing = [name for name in row_of if name not in reached]
        assert missing, "the shape under test covers every reference"

        for field, values in {**per_base(handle), **per_reference(handle)}.items():
            for name in missing:
                assert np.isnan(values[row_of[name]]).all(), \
                    f"{field}: {name} was never named but is not NaN"


def test_a_reference_whose_reads_were_all_turned_away_is_zero(datasets, tmp_path):
    """Nothing counted is not the same as nothing to count. A reference the
    alignments named keeps a zero, however little survived the filter, because
    the reads really were there to be rejected."""
    data = datasets("ragged")
    output = tmp_path / "rejected.h5"
    summary = run_cmuts(data, output, min_mapq=REJECTS_EVERYTHING)

    assert summary.kept == 0, "the filter under test let something through"

    lengths = {name: len(seq) for name, seq in sequences(data.fasta).items()}
    reached = references_with_reads(data.bam)

    with h5py.File(output, "r") as handle:
        row_of = rows_by_name(handle)

        for field, values in per_base(handle).items():
            for name in reached:
                within = values[row_of[name]][:lengths[name]]
                assert not np.isnan(within).any(), f"{field}: {name} went to NaN"
                assert (within == 0).all(), f"{field}: {name} counted something"

        for field, values in per_reference(handle).items():
            for name in reached:
                assert values[row_of[name]] == 0 or field == "reads_filtered", \
                    f"{field}: {name} is not zero"
