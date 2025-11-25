#!/usr/bin/env bash
# Test iHS on PAR1 with constrained parameters to limit haplotype explosion

set -euo pipefail

DATA_HAPBIN_DIR="/faststorage/project/mit-ihh-pib/data/grch38/hapbin"
TEST_DIR="/faststorage/project/mit-ihh-pib/data/grch38/hapbin/test_par1"
REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"

mkdir -p "$TEST_DIR"

CHR="X"
HAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.hap"
MAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.map"

# Test full PAR1 region: positions 288060 - 2781479 (the actual PAR1 boundary)
START_POS=288060
END_POS=2781479

TEST_HAP="$TEST_DIR/chrX_fullPAR1.hap"
TEST_MAP="$TEST_DIR/chrX_fullPAR1.map"
TEST_OUT="$REPO_ROOT/results/test_par1_constrained.ihs.tsv"

echo "Testing iHS on FULL PAR1 with constrained parameters"
echo "Region: ${START_POS}-${END_POS} (actual PAR1 boundary)"
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
echo "  Range: $START_POS - $END_POS (full PAR1)"
echo ""

if [[ "$HAP_LINES" != "$MAP_LINES" ]]; then
    echo "ERROR: Line count mismatch!"
    exit 1
fi

echo "Running iHS with constrained parameters:"
echo "  --max-extend 500000   (limit EHH extension to 500kb)"
echo "  --cutoff 0.10         (higher cutoff = faster EHH decay)"
echo "  --minmaf 0.01         (lower MAF threshold to include more variants)"
echo ""
echo "Start time: $(date)"
echo ""

time pixi run ihsbin \
    --hap "$TEST_HAP" \
    --map "$TEST_MAP" \
    --out "$TEST_OUT" \
    --max-extend 500000 \
    --cutoff 0.10 \
    --minmaf 0.01

echo ""
echo "End time: $(date)"
echo ""

if [[ -s "$TEST_OUT" ]]; then
    OUT_LINES=$(wc -l < "$TEST_OUT")
    VALID_LINES=$((OUT_LINES - 1))
    echo "✓ iHS completed successfully!"
    echo "  Output: $OUT_LINES lines ($VALID_LINES valid results)"
    echo ""
    if [[ $VALID_LINES -gt 0 ]]; then
        echo "First 10 results:"
        head -11 "$TEST_OUT"
        echo ""
        echo "Last 5 results:"
        tail -5 "$TEST_OUT"
    else
        echo "WARNING: No valid iHS results in output (only header)"
    fi
else
    echo "✗ iHS failed - no output created"
    exit 1
fi
