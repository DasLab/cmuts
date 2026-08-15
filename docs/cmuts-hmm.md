# cmuts-hmm

## Purpose

Computing reactivity rates from alignment files via the pair-HMM.

## Requires

- One or more coordinate-sorted alignment files. SAM, BAM, and CRAM formats are all supported
- The FASTA library

## Model Parameters

The pair HMM is configured by five internal parameters. `--dump-params` writes them in the form `--params` reads.

```sh
cmuts-hmm --dump-params > params.txt
```

You may specify a subset of the parameters to modify only them.

## Error Checking

`cmuts-hmm` verifies each FASTA sequence against the alignment header: its name and its length, and, where the header declares an MD5 checksum, the bases themselves. Any mismatch ends the run early. This behavior is configurable via the `--verify` flag.

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
| `--verify CHECKS` | header fields to verify against the FASTA (name, length, checksum, none; default name,length,checksum) |

### Filtering

| Option | Description |
| --- | --- |
| `-q, --min-mapq N` | discard alignments below this mapping quality (0 to 254; default 20) |
| `--min-length N` | discard reads shorter than this (default: no limit) |
| `--max-length N` | discard reads longer than this (default: no limit) |
| `-s, --strand STRANDS` | keep alignments on these strands (forward, reverse; default forward,reverse) |

### Counting

| Option | Description |
| --- | --- |
| `--band N` | reference positions the marginal may look either side of the CIGAR (default 2) |
| `--min-depth D` | evidence a position needs before its rate is written (default 1) |
| `--params FILE` | read the pair HMM's rates from this file (default: built in) |
| `--substitution-weight W` | what a substitution counts towards the mutation total (0 to 1; default 1) |
| `--deletion-weight W` | what a deletion counts towards the mutation total (0 to 1; default 1) |
| `--insertion-weight W` | what an insertion counts towards the mutation total (0 to 1; default 0) |

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
