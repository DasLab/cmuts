"""cmuts-hmm requires a coordinate-sorted file and must say so.

A reference is finished the moment the reader moves past it, which holds only
when reads arrive grouped by reference.
"""

import subprocess

import pytest

from datasets import DATASETS
from support import counted, header_of, reheadered, try_cmuts

# cmuts-hmm reads the sort order from the header alone, so each test declares
# on what the header it was given holds.


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_coordinate_sorted_file_is_accepted(datasets, falsifiable, tmp_path, name):
    data = datasets(name)

    falsifiable("SO:coordinate" in header_of(data))

    assert try_cmuts(data, tmp_path / "out.h5").returncode == 0


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_name_sorted_file_is_refused(datasets, falsifiable, tmp_path, name):
    data = datasets(name)
    bam = tmp_path / "byname.bam"
    with open(bam, "wb") as handle:
        subprocess.run(["samtools", "sort", "-n", str(data.bam)],
                       check=True, stdout=handle, stderr=subprocess.DEVNULL)

    byname = counted((bam,), data.fasta)

    falsifiable(header_of(byname) != header_of(data))

    attempt = try_cmuts(byname, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not coordinate sorted" in attempt.stderr


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_header_without_an_hd_line_is_refused(datasets, falsifiable, tmp_path, name):
    data = datasets(name)
    stripped = reheadered(data, tmp_path,
                          lambda text: "".join(line + "\n" for line in text.splitlines()
                                               if not line.startswith("@HD")))

    falsifiable(header_of(stripped) != header_of(data))

    attempt = try_cmuts(stripped, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not coordinate sorted" in attempt.stderr


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_sort_order_that_merely_starts_the_same_is_refused(datasets, falsifiable,
                                                             tmp_path, name):
    data = datasets(name)
    altered = reheadered(data, tmp_path,
                         lambda text: text.replace("SO:coordinate", "SO:coordinated"))

    falsifiable(header_of(altered) != header_of(data))

    attempt = try_cmuts(altered, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not coordinate sorted" in attempt.stderr


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_other_tags_on_the_hd_line_do_not_hide_the_sort_order(datasets, falsifiable,
                                                              tmp_path, name):
    data = datasets(name)
    altered = reheadered(
        data, tmp_path,
        lambda text: text.replace("@HD\tVN:1.6\tSO:coordinate",
                                  "@HD\tVN:1.6\tSO:coordinate\tGO:none"),
    )

    falsifiable(header_of(altered) != header_of(data))

    assert try_cmuts(altered, tmp_path / "out.h5").returncode == 0


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_an_hd_line_without_a_sort_order_is_refused(datasets, falsifiable, tmp_path,
                                                    name):
    data = datasets(name)
    altered = reheadered(data, tmp_path,
                         lambda text: text.replace("\tSO:coordinate", ""))

    falsifiable(header_of(altered) != header_of(data))

    attempt = try_cmuts(altered, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not coordinate sorted" in attempt.stderr


@pytest.mark.parametrize("name", sorted(DATASETS))
def test_a_sort_order_on_a_line_that_is_not_hd_is_refused(datasets, falsifiable,
                                                          tmp_path, name):
    """The SAM spec carries the sort order on @HD and nowhere else."""
    data = datasets(name)
    def rewrite(text):
        lines = [line for line in text.splitlines() if not line.startswith("@HD")]
        lines[0] = lines[0].replace("@SQ\t", "@SQ\tSO:coordinate\t", 1)
        return "".join(line + "\n" for line in lines)

    altered = reheadered(data, tmp_path, rewrite)

    falsifiable(header_of(altered) != header_of(data))

    attempt = try_cmuts(altered, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not coordinate sorted" in attempt.stderr
