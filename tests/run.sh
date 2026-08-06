#!/usr/bin/env bash
#
# Checks cmuts's filtering against samtools, on alignments written by
# cmuts-gen.
#
# samtools is the oracle throughout: it decides independently how many reads
# should survive each combination of criteria, and the test asks only whether
# cmuts agrees. Nothing here inspects what the processing step computed, so
# these tests stay valid however that changes.
#
# Author: Hamish M. Blair <hmblair@stanford.edu>

set -uo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
GEN=$ROOT/build/cmuts-gen
CMUTS=$ROOT/build/cmuts
SUMMARIZE=$ROOT/tests/summarize.py
WORK=$ROOT/build/tests

FUZZ_ROUNDS=${FUZZ_ROUNDS:-12}

passed=0
failed=0

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

check() {  # what expected actual
    if [ "$2" = "$3" ]; then
        passed=$((passed + 1))
    else
        failed=$((failed + 1))
        printf '  FAIL  %-46s expected %-9s got %s\n' "$1" "$2" "$3"
    fi
}

note() { printf '\n%s\n' "$1"; }

require() {
    for tool in "$@"; do
        command -v "$tool" >/dev/null || { echo "$tool is required"; exit 2; }
    done
}

# ---------------------------------------------------------------------------
# What samtools says
# ---------------------------------------------------------------------------

# Reads surviving a combination of criteria. Strand is a flag samtools knows;
# length it does not, so the sequence column is measured directly. A bound of
# zero is not applied, matching cmuts.
expected_kept() {  # bam mapq strand min max
    local flags=(-F 4 -q "$2")

    case $3 in
        forward) flags+=(-F 16) ;;
        reverse) flags+=(-f 16) ;;
    esac

    samtools view "${flags[@]}" "$1" | awk -v lo="$4" -v hi="$5" '
        { n = length($10)
          if ((lo == 0 || n >= lo) && (hi == 0 || n <= hi)) kept++ }
        END { print kept + 0 }'
}

# Totals that do not depend on the filter, so they are measured once per
# fixture rather than once per case. Held in plain variables, since associative
# arrays would need a bash newer than macOS ships.
MAPPED=0
UNMAPPED=0
TOUCHED=0

measure_fixture() {  # name
    local bam=$WORK/$1.bam

    MAPPED=$(samtools view -c -F 4 "$bam")
    UNMAPPED=$(samtools view -c -f 4 "$bam")
    TOUCHED=$(samtools view -F 4 "$bam" | cut -f3 | sort -u | wc -l | tr -d ' ')
}

# ---------------------------------------------------------------------------
# One combination of criteria
# ---------------------------------------------------------------------------

run_case() {  # fixture mapq strand min max
    local fixture=$1 mapq=$2 strand=$3 lo=$4 hi=$5
    local label="$fixture q=$mapq s=$strand len=$lo-$hi"

    if ! $CMUTS -f "$WORK/$fixture.fasta" -o "$WORK/out.h5" -j 4 \
                -q "$mapq" -s "$strand" --min-length "$lo" --max-length "$hi" \
                "$WORK/$fixture.bam"; then
        failed=$((failed + 1))
        printf '  FAIL  %-46s cmuts exited non-zero\n' "$label"
        return
    fi

    local kept rejected rows unmapped
    read -r kept rejected rows unmapped <<<"$(python3 "$SUMMARIZE" "$WORK/out.h5")"

    check "$label kept" "$(expected_kept "$WORK/$fixture.bam" "$mapq" "$strand" "$lo" "$hi")" "$kept"

    # Every mapped read is either kept or rejected, and a reference that
    # received any mapped read gets a row whether or not anything survived.
    check "$label accounted" "$MAPPED" "$((kept + rejected))"
    check "$label unmapped"  "$UNMAPPED" "$unmapped"
    check "$label rows"      "$TOUCHED" "$rows"
}

# The result must not depend on how many threads produced it.
check_worker_invariance() {  # fixture
    local fixture=$1

    $CMUTS -f "$WORK/$fixture.fasta" -o "$WORK/w1.h5"  -j 1  -q 10 "$WORK/$fixture.bam"
    $CMUTS -f "$WORK/$fixture.fasta" -o "$WORK/w16.h5" -j 16 -q 10 "$WORK/$fixture.bam"

    check "$fixture identical at 1 and 16 workers" "same" \
          "$(python3 - "$WORK/w1.h5" "$WORK/w16.h5" <<'PY'
import sys, h5py, numpy as np
a, b = (h5py.File(p) for p in sys.argv[1:3])
keys = ['coverage', 'mutations', 'reads', 'reads_filtered']
print('same' if all(np.array_equal(a[k][:], b[k][:], equal_nan=True) for k in keys) else 'differ')
PY
)"
}

# The generator is checked too: samtools recomputes MD and NM from the
# alignment and the reference, so a fixture that disagrees is a bad fixture
# rather than a failing filter.
check_fixture_valid() {  # fixture
    local bam=$WORK/$1.bam

    samtools quickcheck "$bam" || check "$1 is well formed" "yes" "no"
    samtools faidx "$WORK/$1.fasta"
    samtools calmd -b "$bam" "$WORK/$1.fasta" >"$WORK/calmd.bam" 2>/dev/null

    local tags='{ for (i = 12; i <= NF; i++) if ($i ~ /^(MD|NM):/) printf "%s ", $i; print "" }'

    if diff -q <(samtools view "$bam" | awk "$tags") \
               <(samtools view "$WORK/calmd.bam" | awk "$tags") >/dev/null; then
        passed=$((passed + 1))
    else
        failed=$((failed + 1))
        printf '  FAIL  %-46s MD/NM disagree with samtools calmd\n' "$1"
    fi
}

# ---------------------------------------------------------------------------
# Fixtures, each chosen for a shape that is easy to get wrong
# ---------------------------------------------------------------------------

make_fixtures() {
    $GEN -o "$WORK/plain"   --seed 101 --references 40  --ref-length 150:600 --reads-per-ref 25
    $GEN -o "$WORK/sparse"  --seed 102 --references 800 --covered 0.3 --reads-per-ref 1:6 \
                            --ref-length 200:300
    $GEN -o "$WORK/single"  --seed 103 --references 1   --ref-length 600 --reads-per-ref 3000
    $GEN -o "$WORK/ragged"  --seed 104 --references 30  --ref-length 60:4000 --reads-per-ref 10:40
    $GEN -o "$WORK/clipped" --seed 105 --references 25  --reads-per-ref 20 \
                            --soft-clips 2 --soft-clip-length 20:80
    $GEN -o "$WORK/indels"  --seed 106 --references 25  --reads-per-ref 20 \
                            --insertions 1:3 --insertion-length 20:200 --deletions 1:2
    $GEN -o "$WORK/lowqual" --seed 107 --references 15  --reads-per-ref 20 --mapq 0
    $GEN -o "$WORK/clean"   --seed 108 --references 15  --reads-per-ref 20 \
                            --mismatch-rate 0 --insertions 0 --deletions 0 --soft-clips 0
}

FIXTURES=(plain sparse single ragged clipped indels lowqual clean)

# ---------------------------------------------------------------------------
# The matrix
# ---------------------------------------------------------------------------

run_matrix() {
    local fixture

    for fixture in "${FIXTURES[@]}"; do
        measure_fixture "$fixture"

        run_case "$fixture"  0 both     0    0      # nothing filtered at all
        run_case "$fixture"  1 both     0    0      # removes the zero-quality reads
        run_case "$fixture" 30 both     0    0
        run_case "$fixture" 60 both     0    0      # above most of the spread
        run_case "$fixture"  0 forward  0    0
        run_case "$fixture"  0 reverse  0    0
        run_case "$fixture"  0 both   200    0      # lower bound only
        run_case "$fixture"  0 both     0  300      # upper bound only
        run_case "$fixture"  0 both   150  400      # a band
        run_case "$fixture"  0 both  9000    0      # admits nothing
        run_case "$fixture" 30 reverse 100 500      # every criterion at once
    done
}

# ---------------------------------------------------------------------------
# Fuzzing, where the case neither of us thought to write turns up
# ---------------------------------------------------------------------------

fuzz() {
    local round seed refs reads mapq strand lo hi

    for ((round = 1; round <= FUZZ_ROUNDS; round++)); do
        seed=$((RANDOM * 32768 + RANDOM))

        refs=$((RANDOM % 60 + 1))
        reads=$((RANDOM % 30 + 1))

        $GEN -o "$WORK/fuzz" --seed "$seed" \
             --references "$refs" \
             --ref-length "$((RANDOM % 200 + 50)):$((RANDOM % 2000 + 300))" \
             --covered "0.$((RANDOM % 9 + 1))" \
             --reads-per-ref "0:$reads" \
             --read-length "$((RANDOM % 80 + 20)):$((RANDOM % 400 + 100))" \
             --mismatch-rate "0.0$((RANDOM % 9))" \
             --insertions "0:$((RANDOM % 3))" --insertion-length "1:$((RANDOM % 150 + 1))" \
             --deletions "0:$((RANDOM % 3))" --deletion-length "1:$((RANDOM % 20 + 1))" \
             --soft-clips "0:2" --soft-clip-length "1:$((RANDOM % 60 + 1))" \
             --unmapped "0:$((RANDOM % 20))" || { echo "  generator failed at seed $seed"; continue; }

        measure_fixture fuzz

        mapq=$((RANDOM % 62))
        lo=$((RANDOM % 400))
        hi=$((RANDOM % 600 + lo))
        case $((RANDOM % 3)) in
            0) strand=both ;;
            1) strand=forward ;;
            *) strand=reverse ;;
        esac

        run_case fuzz "$mapq" "$strand" "$lo" "$hi"
    done
}

# ---------------------------------------------------------------------------

main() {
    require samtools python3
    python3 -c 'import h5py' 2>/dev/null || { echo "h5py is required to read the output"; exit 2; }
    [ -x "$GEN" ] && [ -x "$CMUTS" ] || { echo "run make first"; exit 2; }

    rm -rf "$WORK"
    mkdir -p "$WORK"

    note "Building fixtures"
    make_fixtures

    note "Fixtures agree with samtools calmd"
    for fixture in "${FIXTURES[@]}"; do
        measure_fixture "$fixture"
        check_fixture_valid "$fixture"
    done

    note "Filtering matches samtools"
    run_matrix

    note "Results do not depend on the worker count"
    for fixture in plain sparse single; do
        check_worker_invariance "$fixture"
    done

    note "Fuzzing, $FUZZ_ROUNDS rounds"
    fuzz

    note "$passed passed, $failed failed"
    [ "$failed" -eq 0 ]
}

main "$@"
