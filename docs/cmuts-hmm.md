# cmuts-hmm

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
| `--substitution-weight W` | what a substitution counts towards the mutation total (0 to 1; default 1) |
| `--deletion-weight W` | what a deletion counts towards the mutation total (0 to 1; default 1) |
| `--insertion-weight W` | what an insertion counts towards the mutation total (0 to 1; default 0) |

### Performance

| Option | Description |
| --- | --- |
| `-j, --workers N` | threads running the processing step (default 4) |
| `--decode-threads N` | htslib threads for BGZF decompression (default 4) |
| `--queue-capacity N` | reads in transit at once (default 4096) |
| `--batch N` | reads transferred per queue operation (default 64) |
| `--live-refs N` | references in flight (default 64) |

### Information

| Option | Description |
| --- | --- |
| `-h, --help` | show this help and exit |
| `-V, --version` | show the version and exit |
<!-- END GENERATED cmuts-hmm OPTIONS -->
