#!/usr/bin/env bash
# Test selscan iHS on PAR1 region

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
TEST_DIR="/faststorage/project/mit-ihh-pib/data/grch38/hapbin/test_par1"

cd "$REPO_ROOT"

TEST_HAP="$TEST_DIR/chrX_fullPAR1.hap"
TEST_MAP="$TEST_DIR/chrX_fullPAR1.map"
TEST_OUT="$REPO_ROOT/results/test_par1_selscan"

echo "=========================================="
echo "Testing selscan iHS on full PAR1"
echo "=========================================="
echo "Region: PAR1 (positions 288060 - 2781479)"
echo ""
echo "Input files:"
echo "  HAP: $TEST_HAP"
echo "  MAP: $TEST_MAP"
echo ""

HAP_LINES=$(wc -l < "$TEST_HAP")
echo "  SNPs: $HAP_LINES"
echo ""

echo "Running selscan with parameters:"
echo "  --max-extend 500000  (limit to 500kb)"
echo "  --maf 0.0            (no MAF filter)"
echo "  --threads 8          (parallel processing)"
echo ""
echo "Start time: $(date)"
echo ""

time pixi run selscan \
    --ihs \
    --hap "$TEST_HAP" \
    --map "$TEST_MAP" \
    --out "$TEST_OUT" \
    --max-extend 500000 \
    --maf 0.0 \
    --threads 8

echo ""
echo "End time: $(date)"
echo "=========================================="

if [[ -f "${TEST_OUT}.ihs.out" ]]; then
    OUT_LINES=$(wc -l < "${TEST_OUT}.ihs.out")
    echo "✓ selscan completed!"
    echo "  Output file: ${TEST_OUT}.ihs.out"
    echo "  Lines: $OUT_LINES"
    echo ""

    if [[ $OUT_LINES -gt 1 ]]; then
        echo "First 10 results:"
        head -11 "${TEST_OUT}.ihs.out"
        echo ""
        echo "Last 5 results:"
        tail -5 "${TEST_OUT}.ihs.out"
    fi
else
    echo "✗ No output file created"
    exit 1
fi
