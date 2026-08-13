# Running cmuts

Every run takes the reference sequences and one or more coordinate-sorted alignment files, and writes an HDF5 file. `samtools sort` produces the alignments; no index is needed, and BAM, SAM and CRAM are all read.

## One sample

```sh
cmuts -f references.fasta -o counts.h5 alignments.bam
```

## Reads split across several files

Files are read as one, so a run spread over lanes needs no merging first.

```sh
cmuts -f references.fasta -o counts.h5 lane1.bam lane2.bam lane3.bam
```

## A treated sample and its background

Counting a treated sample and an untreated one gives two files, and `cmuts-sub` takes the background off the signal. The result holds the same datasets, so whatever reads a cmuts output reads it too.

```sh
cmuts -f references.fasta -o treated.h5 treated.bam
cmuts -f references.fasta -o untreated.h5 untreated.bam
cmuts-sub -o reactivity.h5 treated.h5 untreated.h5
```

Add `--clip` to raise a negative reactivity to zero, for the tools downstream that require it.

## Normalizing against a denatured control

A denatured sample measures how reachable each position is with no structure to hide it. Dividing by it puts different references on a comparable scale.

```sh
cmuts-sub -o reactivity.h5 -d denatured.h5 treated.h5 untreated.h5
```

## Choosing which reads count

```sh
cmuts -f references.fasta -o counts.h5 \
      --min-mapq 30 --min-length 100 --max-length 500 --strand forward \
      alignments.bam
```

The length bounds are on the sequence a read stores, which is not the span it covers on the reference: soft-clipped ends and insertions store more than they align to, and deletions store less. Some reads are refused whatever is set here, which [How cmuts counts](counting.md) lists.

## Replacing a result

An existing output file is never overwritten unless `--overwrite` is given, a run costing far more than the command that starts it.
