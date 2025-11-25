#!/usr/bin/env bash
# Test iHS on PAR1 with VERY aggressive parameter constraints

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
TEST_DIR="/faststorage/project/mit-ihh-pib/data/grch38/hapbin/test_par1"

cd "$REPO_ROOT"

TEST_HAP="$TEST_DIR/chrX_fullPAR1.hap"
TEST_MAP="$TEST_DIR/chrX_fullPAR1.map"
TEST_OUT="$REPO_ROOT/results/test_par1_veryaggressive.ihs.tsv"

echo "Testing iHS on full PAR1 with VERY aggressive constraints:"
echo "  --max-extend 100000   (limit to 100kb)"
echo "  --cutoff 0.25         (very high cutoff)"
echo "  --minmaf 0.0          (no MAF filter)"
echo ""
echo "Start time: $(date)"

time pixi run ihsbin \
    --hap "$TEST_HAP" \
    --map "$TEST_MAP" \
    --out "$TEST_OUT" \
    --max-extend 100000 \
    --cutoff 0.25 \
    --minmaf 0.0

echo "End time: $(date)"

if [[ -s "$TEST_OUT" ]]; then
    OUT_LINES=$(wc -l < "$TEST_OUT")
    VALID_LINES=$((OUT_LINES - 1))
    echo "✓ Output: $VALID_LINES valid results"
    if [[ $VALID_LINES -gt 0 ]]; then
        head -6 "$TEST_OUT"
    fi
else
    echo "✗ No output"
fi
