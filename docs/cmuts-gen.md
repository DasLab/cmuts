# cmuts-gen

## CLI Options

<!-- BEGIN GENERATED cmuts-gen OPTIONS -->
### Output

| Option | Description |
| --- | --- |
| `-o, --output PREFIX` | write PREFIX.bam and PREFIX.fasta (required) |
| `--format FORMAT` | alignment format to write (bam, sam; default bam) |

### Layout

| Option | Description |
| --- | --- |
| `--references N` | how many references to write (default 100) |
| `--ref-length DISTRIBUTION` | length of each reference (default 400) |
| `--covered F` | fraction of references receiving any reads (0 to 1; default 1) |
| `--reads-per-ref DISTRIBUTION` | reads on each covered reference (default 20) |

### Reads

| Option | Description |
| --- | --- |
| `--read-length DISTRIBUTION` | reference span each read covers (default 100:400) |
| `--mapq DISTRIBUTION` | mapping quality of each read (default 0,1,10,30,42,60) |
| `--base-quality DISTRIBUTION` | PHRED score of each base (default 30:40) |
| `--reverse F` | fraction of reads on the reverse strand (0 to 1; default 0.5) |
| `--unmapped DISTRIBUTION` | reads aligning nowhere (default 10) |

### Differences from the reference

| Option | Description |
| --- | --- |
| `--mismatch-rate F` | per aligned base (0 to 1; default 0.01) |
| `--insertions DISTRIBUTION` | insertion events per read (default 0:1) |
| `--insertion-length DISTRIBUTION` | bases per insertion (default 1:5) |
| `--deletions DISTRIBUTION` | deletion events per read (default 0:1) |
| `--deletion-length DISTRIBUTION` | bases per deletion (default 1:3) |
| `--soft-clips DISTRIBUTION` | clipped ends per read, none through both (default 0:2) |
| `--soft-clip-length DISTRIBUTION` | bases per clipped end (default 5:30) |

### Determinism

| Option | Description |
| --- | --- |
| `--seed N` | everything generated follows from this (default 1) |

### Information

| Option | Description |
| --- | --- |
| `-h, --help` | show this help and exit |
| `-V, --version` | show the version and exit |

### Advanced

Accepted, and left out of `--help`.

| Option | Description |
| --- | --- |
| `--dump-options` | describe every argument as JSON and exit |
<!-- END GENERATED cmuts-gen OPTIONS -->
