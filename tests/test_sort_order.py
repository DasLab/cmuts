"""cmuts requires a coordinate-sorted file and must say so.

The whole design rests on it: a reference is finished the moment the reader
moves past it, which is only true when reads arrive grouped by reference.
"""

import subprocess

import pytest

from support import try_cmuts, Dataset


def reheadered(data, tmp_path, transform):
    """The same alignments behind a header the transform has rewritten."""
    header = subprocess.run(["samtools", "view", "-H", str(data.bam)],
                            check=True, capture_output=True, text=True).stdout

    path = tmp_path / "header.sam"
    path.write_text(transform(header))

    bam = tmp_path / "reheadered.bam"
    with open(bam, "wb") as handle:
        subprocess.run(["samtools", "reheader", str(path), str(data.bam)],
                       check=True, stdout=handle, stderr=subprocess.DEVNULL)

    return Dataset(bam=bam, fasta=data.fasta, mapped=data.mapped,
                   unmapped=data.unmapped, touched=data.touched)


@pytest.fixture
def data(datasets):
    return datasets("plain")


def test_a_coordinate_sorted_file_is_accepted(data, tmp_path):
    assert try_cmuts(data, tmp_path / "out.h5").returncode == 0


def test_a_name_sorted_file_is_refused(data, tmp_path):
    bam = tmp_path / "byname.bam"
    with open(bam, "wb") as handle:
        subprocess.run(["samtools", "sort", "-n", str(data.bam)],
                       check=True, stdout=handle, stderr=subprocess.DEVNULL)

    attempt = try_cmuts(
        Dataset(bam=bam, fasta=data.fasta, mapped=0, unmapped=0, touched=0),
        tmp_path / "out.h5",
    )

    assert attempt.returncode != 0
    assert "not coordinate sorted" in attempt.stderr


def test_a_header_without_an_hd_line_is_refused(data, tmp_path):
    stripped = reheadered(data, tmp_path,
                          lambda text: "".join(line + "\n" for line in text.splitlines()
                                               if not line.startswith("@HD")))

    assert try_cmuts(stripped, tmp_path / "out.h5").returncode != 0


def test_a_sort_order_that_merely_starts_the_same_is_refused(data, tmp_path):
    """The value runs to the end of its field, so "coordinated" is not it."""
    altered = reheadered(data, tmp_path,
                         lambda text: text.replace("SO:coordinate", "SO:coordinated"))

    assert try_cmuts(altered, tmp_path / "out.h5").returncode != 0


def test_other_tags_on_the_hd_line_do_not_confuse_it(data, tmp_path):
    """SO may sit anywhere among the tags, and other tags may contain its name."""
    altered = reheadered(
        data, tmp_path,
        lambda text: text.replace("@HD\tVN:1.6\tSO:coordinate",
                                  "@HD\tVN:1.6\tSO:coordinate\tGO:none"),
    )

    assert try_cmuts(altered, tmp_path / "out.h5").returncode == 0


def test_an_hd_line_without_a_sort_order_is_refused(data, tmp_path):
    """A header may say nothing about order, which is not the same as saying
    it is sorted."""
    altered = reheadered(data, tmp_path,
                         lambda text: text.replace("\tSO:coordinate", ""))

    assert try_cmuts(altered, tmp_path / "out.h5").returncode != 0


def test_a_sort_order_on_a_line_that_is_not_hd_is_not_believed(data, tmp_path):
    """Only @HD carries a sort order. With no @HD at all, whatever ends up
    first says nothing however it is spelt, which is what the check for the
    line being @HD is there for."""
    def rewrite(text):
        lines = [line for line in text.splitlines() if not line.startswith("@HD")]
        lines[0] = lines[0].replace("@SQ\t", "@SQ\tSO:coordinate\t", 1)
        return "".join(line + "\n" for line in lines)

    assert try_cmuts(reheadered(data, tmp_path, rewrite),
                     tmp_path / "out.h5").returncode != 0
