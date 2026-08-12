"""cmuts requires a coordinate-sorted file and must say so.

A reference is finished the moment the reader moves past it, which holds only
when reads arrive grouped by reference.
"""

import subprocess

import pytest

from support import counted, reheadered, try_cmuts


def test_a_coordinate_sorted_file_is_accepted(data, tmp_path):
    assert try_cmuts(data, tmp_path / "out.h5").returncode == 0


def test_a_name_sorted_file_is_refused(data, tmp_path):
    bam = tmp_path / "byname.bam"
    with open(bam, "wb") as handle:
        subprocess.run(["samtools", "sort", "-n", str(data.bam)],
                       check=True, stdout=handle, stderr=subprocess.DEVNULL)

    attempt = try_cmuts(counted((bam,), data.fasta), tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not coordinate sorted" in attempt.stderr


def test_a_header_without_an_hd_line_is_refused(data, tmp_path):
    stripped = reheadered(data, tmp_path,
                          lambda text: "".join(line + "\n" for line in text.splitlines()
                                               if not line.startswith("@HD")))

    attempt = try_cmuts(stripped, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not coordinate sorted" in attempt.stderr


def test_a_sort_order_that_merely_starts_the_same_is_refused(data, tmp_path):
    altered = reheadered(data, tmp_path,
                         lambda text: text.replace("SO:coordinate", "SO:coordinated"))

    attempt = try_cmuts(altered, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not coordinate sorted" in attempt.stderr


def test_other_tags_on_the_hd_line_do_not_hide_the_sort_order(data, tmp_path):
    altered = reheadered(
        data, tmp_path,
        lambda text: text.replace("@HD\tVN:1.6\tSO:coordinate",
                                  "@HD\tVN:1.6\tSO:coordinate\tGO:none"),
    )

    assert try_cmuts(altered, tmp_path / "out.h5").returncode == 0


def test_an_hd_line_without_a_sort_order_is_refused(data, tmp_path):
    altered = reheadered(data, tmp_path,
                         lambda text: text.replace("\tSO:coordinate", ""))
    attempt = try_cmuts(altered, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not coordinate sorted" in attempt.stderr


def test_a_sort_order_on_a_line_that_is_not_hd_is_refused(data, tmp_path):
    """The SAM spec carries the sort order on @HD and nowhere else."""
    def rewrite(text):
        lines = [line for line in text.splitlines() if not line.startswith("@HD")]
        lines[0] = lines[0].replace("@SQ\t", "@SQ\tSO:coordinate\t", 1)
        return "".join(line + "\n" for line in lines)

    attempt = try_cmuts(reheadered(data, tmp_path, rewrite), tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert "not coordinate sorted" in attempt.stderr
