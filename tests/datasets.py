"""The alignments every test is run against, as parameters for cmuts gen.

Each dataset is a case that is easy to get wrong, and none is written with any
one test in mind. The comment on each describes what makes it awkward.
"""

DATASETS = {
    # The ordinary case, which every contract should hold for.
    "plain":   dict(seed=101, references=40, ref_length="150:600", reads_per_ref=25),

    # Many references, most of them barely covered and some not at all.
    "sparse":  dict(seed=102, references=800, covered=0.3, reads_per_ref="1:6",
                    ref_length="200:300"),

    # Every read on one reference, so every worker writes to one accumulator.
    "single":  dict(seed=103, references=1, ref_length=600, reads_per_ref=3000),

    # Lengths ranging sixtyfold, so most rows are mostly padding.
    "ragged":  dict(seed=104, references=30, ref_length="60:4000", reads_per_ref="10:40"),

    # Reads storing far more than the span they align to, through clipped ends.
    "clipped": dict(seed=105, references=25, reads_per_ref=20, soft_clips=2,
                    soft_clip_length="20:80"),

    # The same, through long insertions, with deletions as well.
    "indels":  dict(seed=106, references=25, reads_per_ref=20, insertions="1:3",
                    insertion_length="20:200", deletions="1:2"),

    # Every read below any useful threshold, so no filter accepts any of them.
    "lowqual": dict(seed=107, references=15, reads_per_ref=20, mapq=0),

    # A mapping quality of 255 alongside one that passes and one that does not.
    # 255 is the SAM spec's "unavailable", which cmuts hmm rejects at every
    # threshold and samtools accepts at every threshold.
    "unavailable": dict(seed=112, references=15, reads_per_ref=20,
                        mapq="0,60,255"),

    # No differences from the reference anywhere, so no mutation is counted.
    "clean":   dict(seed=108, references=15, reads_per_ref=20, mismatch_rate=0,
                    insertions=0, deletions=0, soft_clips=0),

    # Reads running past twice the length of the reference they sit on, which
    # is further than the histogram of lengths reaches.
    "overflowing": dict(seed=109, references=12, ref_length=60, reads_per_ref=20,
                        read_length="40:60", soft_clips="0:2",
                        soft_clip_length="20:150"),

    # Ragged and sparsely covered at once, so a reference with padding and no
    # reads at all is among them.
    "patchy":  dict(seed=110, references=60, ref_length="60:900", covered=0.4,
                    reads_per_ref="5:20"),

    # References of one length, some of them uncovered, so a row holding no
    # reads holds no padding either.
    "flat":    dict(seed=111, references=120, ref_length=300, covered=0.4,
                    reads_per_ref="5:20"),

    # References of a handful of bases, shorter than the reads placed on them.
    "tiny":    dict(seed=113, references=20, ref_length="4:30", reads_per_ref=10,
                    read_length="1:10"),

    # One read per reference, the least depth a rate can be computed from.
    "lonely":  dict(seed=114, references=30, reads_per_ref=1),

    # Every mapped read on the reverse strand, so filtering on one strand accepts all of them
    # and the other none.
    "reversed": dict(seed=115, references=15, reads_per_ref=20, reverse=1.0),

    # An insertion longer than the reference, and a deletion spanning most of
    # the read containing it.
    "gaping":  dict(seed=116, references=12, ref_length=200, reads_per_ref=15,
                    read_length="100:150", deletions=1, deletion_length="50:90",
                    insertions=1, insertion_length="200:400"),

    # One short reference under many reads, so every position is covered
    # hundreds deep.
    "deep":    dict(seed=117, references=1, ref_length=30, reads_per_ref=1500,
                    read_length="5:30"),

    # No reads at all, mapped or not, on references the header still declares.
    "empty":   dict(seed=118, references=10, ref_length=200, covered=0, unmapped=0),
}
