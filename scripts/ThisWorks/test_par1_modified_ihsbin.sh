#!/usr/bin/env bash
# Test modified ihsbin (with larger buffers) on PAR1 region

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
TEST_DIR="/faststorage/project/mit-ihh-pib/data/grch38/hapbin/test_par1"
MODIFIED_IHSBIN="$REPO_ROOT/hapbin_source/build_custom/ihsbin"

cd "$REPO_ROOT"

TEST_HAP="$TEST_DIR/chrX_fullPAR1.hap"
TEST_MAP="$TEST_DIR/chrX_fullPAR1.map"
TEST_OUT="$REPO_ROOT/results/test_par1_modified_ihsbin.ihs.tsv"

echo "=========================================="
echo "Testing MODIFIED ihsbin on full PAR1"
echo "=========================================="
echo "Modified parameters:"
echo "  - Initial maxBreadth: 2000 → 20000 (10x increase)"
echo "  - Reallocation increment: 100 → 1000 (10x increase)"
echo ""
echo "Region: PAR1 (positions 288060 - 2781479)"
echo ""

HAP_LINES=$(wc -l < "$TEST_HAP")
echo "Input: $HAP_LINES SNPs in PAR1"
echo ""

echo "Running with parameters:"
echo "  --max-extend 500000  (limit to 500kb)"
echo "  --minmaf 0.0         (no MAF filter)"
echo ""
echo "Start time: $(date)"
echo ""

time "$MODIFIED_IHSBIN" \
    --hap "$TEST_HAP" \
    --map "$TEST_MAP" \
    --out "$TEST_OUT" \
    --max-extend 500000 \
    --minmaf 0.0

echo ""
echo "End time: $(date)"
echo "=========================================="

if [[ -s "$TEST_OUT" ]]; then
    OUT_LINES=$(wc -l < "$TEST_OUT")
    VALID_LINES=$((OUT_LINES - 1))
    echo "✓ Modified ihsbin completed!"
    echo "  Output: $OUT_LINES lines ($VALID_LINES valid results)"
    echo ""
    if [[ $VALID_LINES -gt 0 ]]; then
        echo "SUCCESS! We got valid iHS results for PAR1!"
        echo ""
        echo "First 10 results:"
        head -11 "$TEST_OUT"
        echo ""
        echo "Last 5 results:"
        tail -5 "$TEST_OUT"
    else
        echo "Still no valid results, but completed without hanging"
    fi
else
    echo "✗ No output file created"
    exit 1
fi
