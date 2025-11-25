#!/usr/bin/env bash
# Test iHS on a small slice of PAR1 to see if it can complete

set -euo pipefail

DATA_HAPBIN_DIR="/faststorage/project/mit-ihh-pib/data/grch38/hapbin"
TEST_DIR="/faststorage/project/mit-ihh-pib/data/grch38/hapbin/test_par1"
REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"

mkdir -p "$TEST_DIR"

CHR="X"
HAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.hap"
MAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.map"

# Test a very small slice of PAR1: first 100kb (positions 288060 - 388060)
START_POS=288060
END_POS=388060

TEST_HAP="$TEST_DIR/chrX_par1_test.hap"
TEST_MAP="$TEST_DIR/chrX_par1_test.map"
TEST_OUT="$REPO_ROOT/results/test_par1_small.ihs.tsv"

echo "Testing iHS on small PAR1 slice: ${START_POS}-${END_POS}"
echo ""

# Extract MAP lines for this range
awk -v start="$START_POS" -v end="$END_POS" '$4 >= start && $4 <= end' "$MAP" > "$TEST_MAP"

# Get line numbers
FIRST_LINE=$(awk -v start="$START_POS" '$4 >= start {print NR; exit}' "$MAP")
LAST_LINE=$(awk -v end="$END_POS" '$4 <= end {line=NR} END {print line}' "$MAP")

echo "Extracting lines ${FIRST_LINE}-${LAST_LINE} from HAP file..."
sed -n "${FIRST_LINE},${LAST_LINE}p" "$HAP" > "$TEST_HAP"

# Verify
HAP_LINES=$(wc -l < "$TEST_HAP")
MAP_LINES=$(wc -l < "$TEST_MAP")
HAP_SIZE=$(du -h "$TEST_HAP" | cut -f1)

echo "Test data created:"
echo "  HAP: $HAP_SIZE, $HAP_LINES SNPs"
echo "  MAP: $MAP_LINES positions"
echo "  Range: $START_POS - $END_POS (100kb slice of PAR1)"
echo ""

if [[ "$HAP_LINES" != "$MAP_LINES" ]]; then
    echo "ERROR: Line count mismatch!"
    exit 1
fi

echo "Running iHS on test slice..."
echo "Start time: $(date)"
echo ""

time pixi run ihsbin --hap "$TEST_HAP" --map "$TEST_MAP" --out "$TEST_OUT"

echo ""
echo "End time: $(date)"
echo ""

if [[ -s "$TEST_OUT" ]]; then
    OUT_LINES=$(wc -l < "$TEST_OUT")
    echo "✓ iHS completed successfully!"
    echo "  Output: $OUT_LINES lines"
    echo ""
    echo "First 10 results:"
    head -11 "$TEST_OUT"
else
    echo "✗ iHS failed - no output created"
    exit 1
fi
