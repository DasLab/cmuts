"""The alignments the tests are run against, as parameters for the generator.

Named here rather than in conftest.py, which pytest collects rather than
imports, so that a test may take the whole catalogue by importing it.
"""

# Each dataset is one that is easy to get wrong: many references barely covered,
# one reference deeply covered, lengths ranging sixtyfold, reads whose stored
# length far exceeds the span they align to through soft clipping and again
# through long insertions, a file where every read falls below any useful
# threshold, one with no differences from the reference at all, one whose reads
# run past twice the length of the reference they are placed on, one carrying
# the mapping quality that stands for none at all, and one that is both ragged
# and sparsely covered, so that a reference with padding and no reads at all is
# among them.
DATASETS = {
    "plain":   dict(seed=101, references=40, ref_length="150:600", reads_per_ref=25),
    "sparse":  dict(seed=102, references=800, covered=0.3, reads_per_ref="1:6",
                    ref_length="200:300"),
    "single":  dict(seed=103, references=1, ref_length=600, reads_per_ref=3000),
    "ragged":  dict(seed=104, references=30, ref_length="60:4000", reads_per_ref="10:40"),
    "clipped": dict(seed=105, references=25, reads_per_ref=20, soft_clips=2,
                    soft_clip_length="20:80"),
    "indels":  dict(seed=106, references=25, reads_per_ref=20, insertions="1:3",
                    insertion_length="20:200", deletions="1:2"),
    "lowqual": dict(seed=107, references=15, reads_per_ref=20, mapq=0),
    # 255 alongside qualities that pass and qualities that do not, so that a
    # run over it counts something and refuses something.
    "unavailable": dict(seed=112, references=15, reads_per_ref=20,
                        mapq="0,60,255"),
    "clean":   dict(seed=108, references=15, reads_per_ref=20, mismatch_rate=0,
                    insertions=0, deletions=0, soft_clips=0),
    "overflowing": dict(seed=109, references=12, ref_length=60, reads_per_ref=20,
                        read_length="40:60", soft_clips="0:2",
                        soft_clip_length="20:150"),
    "patchy":  dict(seed=110, references=60, ref_length="60:900", covered=0.4,
                    reads_per_ref="5:20"),
    "flat":    dict(seed=111, references=120, ref_length=300, covered=0.4,
                    reads_per_ref="5:20"),
}
