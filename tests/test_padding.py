"""What a value of NaN means in the output.

The arrays with a row per reference are rectangular and the references are not,
so every row runs to the width the longest reference needs. NaN is what tells
the difference: a column past what its own reference needs, and a reference no
alignment ever named.

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

from support import references_with_reads, rows_by_name, run_cmuts, sequences

# Higher than any mapping quality a read can carry, so every read is turned
# away and the references themselves are all that is left.
REJECTS_EVERYTHING = 61


# Two kinds of row are written, indexed by different things, so which one an
# array is has to be named rather than read off its shape. Everything not
# listed here is indexed by reference position, and runs only as far as its own
# reference; a length is not a position, so those rows are data throughout.
PER_LENGTH = ("read_lengths",)


def rectangular(output):
    """The arrays with a row per reference, whatever indexes the row."""
    return {name: output[name][:] for name in output if output[name].ndim == 2}


def per_base(output):
    """The arrays a reference's own length bounds, which are the padded ones."""
    return {k: v for k, v in rectangular(output).items() if k not in PER_LENGTH}


def row_extent(field, ref_len, width):
    """Columns of a row that hold data; the rest is padding."""
    return width if field in PER_LENGTH else ref_len


def per_reference(output):
    """The arrays with one value for each reference."""
    return {
        name: output[name][:]
        for name in output
        if output[name].ndim == 1 and np.issubdtype(output[name].dtype, np.floating)
    }


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

        for field, values in per_base(handle).items():
            width = values.shape[1]
            shorter = [n for n, ln in lengths.items() if ln < width]
            assert shorter, f"{field}: no reference is narrower than the widest"

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

        for field, values in rectangular(handle).items():
            width = values.shape[1]
            for name in reached:
                within = values[row_of[name]][:row_extent(field, lengths[name], width)]
                assert not np.isnan(within).any(), f"{field}: {name} has a hole in it"


# ---------------------------------------------------------------------------
# A reference no alignment named
# ---------------------------------------------------------------------------


def test_a_reference_no_read_named_is_zero_over_its_own_bases(datasets, tmp_path):
    """Ragged and sparsely covered, so a reference no alignment named has
    padding of its own to be told apart from what it counted.

    Nothing counted is what a zero says; a column outside the reference is what
    NaN says. A reference no read named is the first over all of its bases and
    the second past them, which is what one whose reads were all rejected also
    holds."""
    data = datasets("patchy")
    output = tmp_path / "patchy.h5"
    run_cmuts(data, output)

    lengths = {name: len(seq) for name, seq in sequences(data.fasta).items()}
    widest = max(lengths.values())
    reached = references_with_reads(data.bam)

    with h5py.File(output, "r") as handle:
        row_of = rows_by_name(handle)
        missing = [name for name in row_of if name not in reached]
        assert missing, "the shape under test covers every reference"
        assert any(lengths[name] < widest for name in missing), \
            "no uncovered reference is short enough to carry padding"

        for field, values in rectangular(handle).items():
            width = values.shape[1]
            for name in missing:
                row = values[row_of[name]]
                extent = row_extent(field, lengths[name], width)

                assert (row[:extent] == 0).all(), \
                    f"{field}: {name} was named by nothing and is not zero"
                assert np.isnan(row[extent:]).all(), \
                    f"{field}: {name} has padding that is not NaN"

        # A count of reads is zero where no read arrived, there being nothing
        # a count could mean by NaN that zero does not say.
        for field, values in per_reference(handle).items():
            for name in missing:
                assert values[row_of[name]] == 0, \
                    f"{field}: {name} was never named but is not zero"


def test_an_uncovered_reference_of_full_length_holds_no_nan(datasets, tmp_path):
    """One length throughout, which is the ordinary shape of a library, so no
    reference has padding and an uncovered row is what the fill alone says."""
    data = datasets("flat")
    output = tmp_path / "flat.h5"
    run_cmuts(data, output)

    reached = references_with_reads(data.bam)

    with h5py.File(output, "r") as handle:
        row_of = rows_by_name(handle)
        missing = [name for name in row_of if name not in reached]
        assert missing, "the shape under test covers every reference"

        for field, values in rectangular(handle).items():
            for name in missing:
                row = values[row_of[name]]
                assert not np.isnan(row).any(), f"{field}: {name} holds a NaN"
                assert (row == 0).all(), f"{field}: {name} is not zero"


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

        for field, values in rectangular(handle).items():
            width = values.shape[1]
            for name in reached:
                within = values[row_of[name]][:row_extent(field, lengths[name], width)]
                assert not np.isnan(within).any(), f"{field}: {name} went to NaN"
                assert (within == 0).all(), f"{field}: {name} counted something"

        for field, values in per_reference(handle).items():
            for name in reached:
                assert values[row_of[name]] == 0 or field == "reads_filtered", \
                    f"{field}: {name} is not zero"
