"""cmuts-hmm requires a coordinate-sorted file and must say so.

A reference is finished the moment the reader moves past it, which holds only
when reads arrive grouped by reference. The sort order is read from the header
alone, so each test declares on what the header it was given holds.
"""

from alignments import measure_dataset, replace_header
from oracle import header_text
from programs import samtools_into, try_cmuts

UNSORTED_MESSAGE = "not coordinate sorted"


def _without_hd_line(text: str) -> str:
    return "".join(line + "\n" for line in text.splitlines()
                   if not line.startswith("@HD"))


def _sort_order_on_the_first_sq_line(text: str) -> str:
    lines = [line for line in text.splitlines() if not line.startswith("@HD")]
    lines[0] = lines[0].replace("@SQ\t", "@SQ\tSO:coordinate\t", 1)

    return "".join(line + "\n" for line in lines)


def test_a_coordinate_sorted_file_is_accepted(data, falsifiable, tmp_path):
    falsifiable("SO:coordinate" in header_text(data.bam))

    assert try_cmuts(data, tmp_path / "out.h5").returncode == 0


def test_a_name_sorted_file_is_refused(data, falsifiable, tmp_path):
    bam = samtools_into(tmp_path / "byname.bam", "sort", "-n", data.bam)
    byname = measure_dataset(data.name, (bam,), data.fasta)

    falsifiable(header_text(byname.bam) != header_text(data.bam))

    attempt = try_cmuts(byname, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert UNSORTED_MESSAGE in attempt.stderr


def test_a_header_without_an_hd_line_is_refused(data, falsifiable, tmp_path):
    stripped = replace_header(data, tmp_path, _without_hd_line)

    falsifiable(header_text(stripped.bam) != header_text(data.bam))

    attempt = try_cmuts(stripped, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert UNSORTED_MESSAGE in attempt.stderr


def test_a_sort_order_that_merely_starts_the_same_is_refused(data, falsifiable,
                                                             tmp_path):
    altered = replace_header(
        data, tmp_path, lambda text: text.replace("SO:coordinate", "SO:coordinated"))

    falsifiable(header_text(altered.bam) != header_text(data.bam))

    attempt = try_cmuts(altered, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert UNSORTED_MESSAGE in attempt.stderr


def test_other_tags_on_the_hd_line_do_not_hide_the_sort_order(data, falsifiable,
                                                              tmp_path):
    altered = replace_header(
        data, tmp_path,
        lambda text: text.replace("@HD\tVN:1.6\tSO:coordinate",
                                  "@HD\tVN:1.6\tSO:coordinate\tGO:none"),
    )

    falsifiable(header_text(altered.bam) != header_text(data.bam))

    assert try_cmuts(altered, tmp_path / "out.h5").returncode == 0


def test_an_hd_line_without_a_sort_order_is_refused(data, falsifiable, tmp_path):
    altered = replace_header(
        data, tmp_path, lambda text: text.replace("\tSO:coordinate", ""))

    falsifiable(header_text(altered.bam) != header_text(data.bam))

    attempt = try_cmuts(altered, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert UNSORTED_MESSAGE in attempt.stderr


def test_a_sort_order_on_a_line_that_is_not_hd_is_refused(data, falsifiable, tmp_path):
    """The SAM spec puts the sort order on @HD and nowhere else."""
    altered = replace_header(data, tmp_path, _sort_order_on_the_first_sq_line)

    falsifiable(header_text(altered.bam) != header_text(data.bam))

    attempt = try_cmuts(altered, tmp_path / "out.h5")

    assert attempt.returncode != 0
    assert UNSORTED_MESSAGE in attempt.stderr
