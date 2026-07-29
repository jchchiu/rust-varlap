#!/usr/bin/env bash
# compare_outputs.sh — compare two CSVs field-by-field, skipping the header row.
# Writes any differences to a diff file and exits non-zero if any are found.
#
# Usage:
#   compare_outputs.sh <file1.csv> <file2.csv> <diff_output.txt>
#
# Exit codes:
#   0 — files match (diff output is empty)
#   1 — differences found (diff output contains details)
#   2 — wrong number of arguments or input file not found

set -euo pipefail

# ── Argument handling ────────────────────────────────────────────────────────

if [[ $# -ne 3 ]]; then
    echo "Usage: $(basename "$0") <file1.csv> <file2.csv> <diff_output.txt>" >&2
    exit 2
fi

FILE1="$1"
FILE2="$2"
DIFF_OUT="$3"

if [[ ! -f "$FILE1" ]]; then
    echo "Error: file1 not found: $FILE1" >&2
    exit 2
fi

if [[ ! -f "$FILE2" ]]; then
    echo "Error: file2 not found: $FILE2" >&2
    exit 2
fi

# ── Comparison (NR==1 skips the header row in both files) ───────────────────

awk -F',' '
NR==FNR {
    if (NR == 1) next
    for (i=1; i<=NF; i++) a[FNR,i]=$i
    max_nf[FNR]=NF
    next
}
FNR == 1 { next }
{
    max = (NF > max_nf[FNR] ? NF : max_nf[FNR])
    for (i=1; i<=max; i++) {
        if (a[FNR,i] != $i) {
            print "Row " FNR ", Col " i ": file1=\"" a[FNR,i] "\", file2=\"" $i "\""
        }
    }
}
' "$FILE1" "$FILE2" > "$DIFF_OUT"

# ── Result ───────────────────────────────────────────────────────────────────

if [[ -s "$DIFF_OUT" ]]; then
    echo "❌ Differences found:" >&2
    cat "$DIFF_OUT" >&2
    exit 1
else
    echo "✅ Outputs match (header row skipped)."
    exit 0
fi