#!/usr/bin/env bash
# Merge iHS results from chromosome X regions into a single file

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
REGION_DIR="$REPO_ROOT/results/ihs_fixed_par1_paraproc/regions"
OUT_DIR="$REPO_ROOT/results/ihs_fixed_par1_paraproc"

FINAL_OUTPUT="$OUT_DIR/ALL.chrX.ihs.tsv"
TEMP_OUTPUT="$OUT_DIR/ALL.chrX.ihs.tmp"

echo "=========================================="
echo "Merging chromosome X iHS results"
echo "=========================================="
echo ""

# Check that all region files exist
NUM_REGIONS=5
MISSING=0

for i in $(seq 1 $NUM_REGIONS); do
    REGION_FILE="$REGION_DIR/chrX_region${i}.ihs.tsv"
    if [[ ! -s "$REGION_FILE" ]]; then
        echo "ERROR: Missing or empty file: $REGION_FILE"
        MISSING=1
    else
        LINES=$(wc -l < "$REGION_FILE")
        SIZE=$(du -h "$REGION_FILE" | cut -f1)
        echo "Region $i: $LINES lines, $SIZE"
    fi
done

echo ""

if [[ $MISSING -eq 1 ]]; then
    echo "ERROR: Some region files are missing. Cannot merge."
    exit 1
fi

echo "Merging files..."

# Extract header from first file
head -1 "$REGION_DIR/chrX_region1.ihs.tsv" > "$TEMP_OUTPUT"

# Append data from all regions (skip headers)
for i in $(seq 1 $NUM_REGIONS); do
    echo "  Adding region $i..."
    tail -n +2 "$REGION_DIR/chrX_region${i}.ihs.tsv" >> "$TEMP_OUTPUT"
done

# Sort by position (column 2)
echo ""
echo "Sorting by genomic position..."
(head -1 "$TEMP_OUTPUT" && tail -n +2 "$TEMP_OUTPUT" | sort -t':' -k2,2n) > "$FINAL_OUTPUT"

# Clean up temp file
rm "$TEMP_OUTPUT"

echo ""
echo "=========================================="
echo "✓ Merge complete"
echo "=========================================="
echo ""

TOTAL_LINES=$(wc -l < "$FINAL_OUTPUT")
TOTAL_SIZE=$(du -h "$FINAL_OUTPUT" | cut -f1)

echo "Final output: $FINAL_OUTPUT"
echo "  Total lines: $TOTAL_LINES (including header)"
echo "  Total variants: $((TOTAL_LINES - 1))"
echo "  File size: $TOTAL_SIZE"
echo ""

# Show first and last few results
echo "First 5 iHS results:"
head -6 "$FINAL_OUTPUT"
echo ""

echo "Last 5 iHS results:"
tail -5 "$FINAL_OUTPUT"
echo ""

# Verify position ordering
echo "Verifying sort order..."
FIRST_POS=$(tail -n +2 "$FINAL_OUTPUT" | head -1 | awk '{print $2}')
LAST_POS=$(tail -1 "$FINAL_OUTPUT" | awk '{print $2}')
echo "  Position range: $FIRST_POS - $LAST_POS"
echo ""

echo "✓ All done!"
