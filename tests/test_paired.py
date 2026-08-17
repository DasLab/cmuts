"""cmuts-hmm requires merged reads and must refuse paired ones.

Two mates read one molecule, so counting them as separate reads would count
their overlap twice. Merging happens before alignment, so a record still
carrying the paired flag marks a file the merge never saw, and the run must
stop rather than count it.
"""

from alignments import mark_paired
from programs import try_cmuts


def test_a_paired_read_is_refused(data, falsifiable, tmp_path):
    """A file holding any paired read is refused, naming the file and the fix.
    A file holding no read holds no pair, and passes."""
    paired = mark_paired(data, tmp_path)
    has_reads = paired.mapped + paired.unmapped > 0

    falsifiable(has_reads)

    attempt = try_cmuts(paired, tmp_path / "out.h5")

    assert (attempt.returncode != 0) == has_reads

    if has_reads:
        assert str(paired.bam) in attempt.stderr
        assert "merge" in attempt.stderr
