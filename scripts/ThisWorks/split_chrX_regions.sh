#!/usr/bin/env bash
# Split chromosome X into genomic regions for parallel iHS computation

set -euo pipefail

DATA_HAPBIN_DIR="/faststorage/project/mit-ihh-pib/data/grch38/hapbin"
OUT_DIR="/faststorage/project/mit-ihh-pib/data/grch38/hapbin/chrX_regions"

mkdir -p "$OUT_DIR"

CHR="X"
HAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.hap"
MAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.map"

# Define 5 genomic regions for chrX (approximate)
# Based on position range: 288060 - 156027371
# Each region ~31M bp

declare -a REGIONS=(
    "1:288060:31500000"
    "2:31500001:62500000"
    "3:62500001:93500000"
    "4:93500001:124500000"
    "5:124500001:156027371"
)

echo "Splitting chromosome X into ${#REGIONS[@]} regions..."
echo ""

for region in "${REGIONS[@]}"; do
    IFS=':' read -r region_num start_pos end_pos <<< "$region"

    OUT_HAP="$OUT_DIR/chrX_region${region_num}.hap"
    OUT_MAP="$OUT_DIR/chrX_region${region_num}.map"

    echo "Region ${region_num}: ${start_pos}-${end_pos}"

    # Extract lines from MAP file within position range
    # MAP format: chr rsid genetic_pos physical_pos
    awk -v start="$start_pos" -v end="$end_pos" '$4 >= start && $4 <= end' "$MAP" > "$OUT_MAP"

    # Get line numbers for this region
    FIRST_LINE=$(awk -v start="$start_pos" '$4 >= start {print NR; exit}' "$MAP")
    LAST_LINE=$(awk -v end="$end_pos" '$4 <= end {line=NR} END {print line}' "$MAP")

    if [[ -z "$FIRST_LINE" || -z "$LAST_LINE" ]]; then
        echo "  WARNING: No SNPs found in this region"
        continue
    fi

    # Extract corresponding lines from HAP file
    sed -n "${FIRST_LINE},${LAST_LINE}p" "$HAP" > "$OUT_HAP"

    # Verify
    HAP_LINES=$(wc -l < "$OUT_HAP")
    MAP_LINES=$(wc -l < "$OUT_MAP")

    echo "  Extracted ${MAP_LINES} SNPs (lines ${FIRST_LINE}-${LAST_LINE})"
    echo "  Output: $(du -h "$OUT_HAP" | cut -f1) HAP, $(du -h "$OUT_MAP" | cut -f1) MAP"

    if [[ "$HAP_LINES" != "$MAP_LINES" ]]; then
        echo "  ERROR: Line count mismatch!"
        exit 1
    fi

    echo ""
done

echo "✓ Chromosome X split into ${#REGIONS[@]} regions"
echo "Output directory: $OUT_DIR"
