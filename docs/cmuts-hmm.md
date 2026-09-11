# `cmuts hmm`

## Purpose

Computing reactivity rates from alignment files via the pair HMM.

## Requires

- One or more coordinate-sorted alignment files of single-end or merged reads. SAM, BAM, and CRAM formats are all supported
- The FASTA library

Alignments must be single-end to avoid double-counting; a read carrying the paired flag causes the program to exit early. Please merge reads upstream before passing data to `cmuts hmm`.

## Model Parameters

The pair HMM is configured by five internal parameters. `--dump-params` writes them in the form `--params` reads.

```sh
cmuts hmm --dump-params > params.txt
```

You may specify a subset of the parameters to modify only them.

## Error Checking

`cmuts hmm` verifies the FASTA against the alignment header, by comparing each sequence's name and length, and its MD5 checksum when present. Any mismatch between the header and the FASTA ends the run early. This behavior is configurable via the `--verify` flag.

```{note}
The length check is required to avoid buffer overflows and cannot be disabled.
```

Both the DNA and the RNA forms of the FASTA sequence are accepted during the checksum comparison.

## Error Computation

`cmuts hmm` computes a per-base estimate for the reactivity error. This is purely the statistical error introduced by finite read depths; it does not account for experimental or systemic errors.

## Rejected Reads

Reads which do not pass the configured filters are rejected, meaning they are not processed by the HMM and do not contribute to the final reactivity. A secondary alignment, a read with MAPQ 255, a read missing a sequence, or a read with no CIGAR is automatically rejected.

A read to which the HMM cannot assign a nonzero alignment probability is also rejected; this can only occur if you set one or more transition probabilities to zero.

Unmapped reads are rejected by nature of having no reference to compute mutation rates against.

## Split Reads

A read the aligner could not place in one piece is written as a primary alignment and one supplementary alignment for each further piece. Every piece is processed and adds to the coverage and the rates, but the totals in the `reads` group count the read once, at its primary alignment, and bin it by the length the primary stores.

Pass `--drop-supplementary` to process the primary alignment alone. This suits a library whose molecules hold more than one copy of the reference, where the pieces cover the same positions and counting them all measures one molecule several times. Leave it unset where reads are longer than the reference and each piece covers positions the others do not.

## Output

<!-- BEGIN GENERATED cmuts-hmm DATASETS -->
| Dataset | Source |
| --- | --- |
| [`mismatches/rate`](format.md#mismatchesrate) | Estimated by the HMM. |
| [`mismatches/error`](format.md#mismatcheserror) | The binomial standard error of the rate at the position's depth. |
| [`insertions/rate`](format.md#insertionsrate) | Estimated by the HMM. |
| [`insertions/error`](format.md#insertionserror) | The binomial standard error of the rate at the position's depth. |
| [`deletions/rate`](format.md#deletionsrate) | Estimated by the HMM. |
| [`deletions/error`](format.md#deletionserror) | The binomial standard error of the rate at the position's depth. |
| [`coverage`](format.md#coverage) | Estimated by the HMM. |
| [`sequence`](format.md#sequence) | Tokenized from the FASTA. |
| [`reads/lengths`](format.md#readslengths) | The read length reported in the alignment. |
| [`reads/counted`](format.md#readscounted) | The number of reads the HMM successfully processed. |
| [`reads/rejected`](format.md#readsrejected) | The number of reads rejected by a filter or by the HMM. |
| [`reads/unmapped`](format.md#readsunmapped) | The number of unmapped reads in the alignment. |
| [`pairwise/correlation`](format.md#pairwisecorrelation) | Estimated by the HMM. |
| [`pairwise/conditional`](format.md#pairwiseconditional) | Estimated by the HMM. |
| [`pairwise/coverage`](format.md#pairwisecoverage) | Estimated by the HMM. |
<!-- END GENERATED cmuts-hmm DATASETS -->

## CLI Options

<!-- BEGIN GENERATED cmuts-hmm OPTIONS -->
### Arguments

| Argument | Description |
| --- | --- |
| `BAM...` | coordinate-sorted alignments |

### Input and output

| Option | Description |
| --- | --- |
| `-f, --fasta FASTA` | reference sequences (required) |
| `-o, --output HDF5` | write results to this file (required) |
| `--overwrite` | replace the output file if it already exists |
| `--verify CHECKS` | identity checks to make against the FASTA (name, checksum, none; default name,checksum) |

### Filtering

| Option | Description |
| --- | --- |
| `-q, --min-mapq N` | reject alignments below this mapping quality (0 to 254; default 20) |
| `--min-length N` | reject reads shorter than this (default: no limit) |
| `--max-length N` | reject reads longer than this (default: no limit) |
| `--drop-supplementary` | reject all but the primary piece of a split read |
| `-s, --strand STRANDS` | reject alignments not on these strands (forward, reverse; default forward,reverse) |

### Counting

| Option | Description |
| --- | --- |
| `--band N` | reference positions the pair HMM may look either side of the CIGAR (default 2) |
| `--min-phred Q` | assign the maximum sequencing error to bases below this quality (0 to 255; default 0) |
| `--pairwise STATS` | write these statistics of how often two positions are modified together (correlation, conditional, none; default none) |
| `--min-depth D` | write NaN rates and errors for positions below this depth (default 1) |
| `--nan-5p N` | write NaN rates and errors for this many bases at the 5' end (default 0) |
| `--nan-3p N` | write NaN rates and errors for this many bases at the 3' end (default 0) |
| `--params FILE` | read the pair HMM's rates from this file (default: built in) |

### Performance

| Option | Description |
| --- | --- |
| `-j, --workers N` | threads running the processing step (default 1) |
| `--decode-threads N` | htslib threads for BGZF decompression (default 0) |

### Information

| Option | Description |
| --- | --- |
| `-h, --help` | show this help and exit |
| `-V, --version` | show the version and exit |
| `--dump-params` | write the rates in the form --params reads and exit |

### Advanced

Accepted, and left out of `--help`.

| Option | Description |
| --- | --- |
| `--queue-capacity N` | reads in transit at once (default 4096) |
| `--batch N` | reads transferred per queue operation (default 64) |
| `--live-refs N` | references in flight (default 64) |
| `--dump-options` | describe every argument as JSON and exit |
| `--dump-layout` | describe the output format as JSON and exit |
<!-- END GENERATED cmuts-hmm OPTIONS -->
