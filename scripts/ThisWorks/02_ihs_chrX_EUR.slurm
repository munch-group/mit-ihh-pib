#!/usr/bin/env bash
#SBATCH -J ihs-chrX-EUR
#SBATCH -A mit-ihh-pib
#SBATCH -p normal
#SBATCH -c 8
#SBATCH -t 36:00:00
#SBATCH -o /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_chrX_EUR_%j.out
#SBATCH -e /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_chrX_EUR_%j.err

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
DATA_DIR="/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_EUR"

cd "$REPO_ROOT"

CHR="X"

OUT_DIR="$REPO_ROOT/results/ihs_EUR"
mkdir -p "$OUT_DIR" "$REPO_ROOT/logs"

HAP="$DATA_DIR/chr${CHR}_EUR.hapbin.hap"
MAP="$DATA_DIR/chr${CHR}_EUR.hapbin.map"
OUT="$OUT_DIR/EUR.chr${CHR}.ihs.tsv"

echo "=========================================="
echo "iHS calculation for chromosome ${CHR}"
echo "Population: EUR ONLY"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================="
echo "Inputs:"
echo "  HAP: $HAP"
echo "  MAP: $MAP"
echo "Output:"
echo "  OUT: $OUT"
echo ""

# Check if output already exists
if [[ -s "$OUT" ]]; then
    echo "✓ Output already exists: $OUT"
    echo "  If you want to rerun, delete this file first."
    exit 0
fi

# Validate inputs
if [[ ! -s "$HAP" ]]; then
    echo "ERROR: HAP file missing or empty: $HAP"
    echo "Did you run 01_subset_EUR_samples_hapbin.py first?"
    exit 1
fi

if [[ ! -s "$MAP" ]]; then
    echo "ERROR: MAP file missing or empty: $MAP"
    exit 1
fi

# Check HAP file format (EUR subset is plain text, not gzipped)
echo "Verifying HAP file format..."
HAP_LINES=$(wc -l < "$HAP")
HAP_SIZE=$(du -h "$HAP" | cut -f1)
FIELDS_LINE1=$(head -1 "$HAP" | awk '{print NF}')
echo "  HAP file: $HAP_SIZE, $HAP_LINES variants"
echo "  First line: $FIELDS_LINE1 haplotypes (EUR samples only)"
echo ""

# Expected haplotypes for EUR:
# 503 samples, but X chromosome varies (males=1, females=2 haplotypes)
# Should be between 503 and 1006 haplotypes
if [[ $FIELDS_LINE1 -lt 500 ]] || [[ $FIELDS_LINE1 -gt 1100 ]]; then
    echo "WARNING: Unexpected number of haplotypes"
    echo "  Expected ~503-1006 for EUR X chromosome"
    echo "  Got: $FIELDS_LINE1"
fi

# Check MAP file
echo "Verifying MAP file..."
MAP_LINES=$(wc -l < "$MAP")
MAP_FIRST_POS=$(head -1 "$MAP" | awk '{print $4}')
MAP_LAST_POS=$(tail -1 "$MAP" | awk '{print $4}')
echo "  MAP file: $MAP_LINES positions"
echo "  Position range: $MAP_FIRST_POS - $MAP_LAST_POS"
echo ""

if [[ "$HAP_LINES" != "$MAP_LINES" ]]; then
    echo "WARNING: HAP and MAP have different line counts"
    echo "  HAP: $HAP_LINES"
    echo "  MAP: $MAP_LINES"
fi

echo "Running iHS for chr${CHR} (EUR samples only)..."
echo "This may take several hours..."
echo ""

pixi run ihsbin --hap "$HAP" --map "$MAP" --out "$OUT"

echo ""
echo "=========================================="
if [[ -s "$OUT" ]]; then
    OUT_LINES=$(wc -l < "$OUT")
    OUT_SIZE=$(du -h "$OUT" | cut -f1)
    echo "✓ iHS completed for chr${CHR} (EUR)"
    echo "=========================================="
    echo ""
    echo "Output file: $OUT"
    echo "  Size: $OUT_SIZE"
    echo "  Lines: $OUT_LINES"
    echo ""

    # Show first few results
    echo "First 5 iHS results:"
    head -6 "$OUT"
    echo ""

    # Statistics
    echo "Quick statistics:"
    tail -n +2 "$OUT" | awk '{
        if ($5 != "NA" && $5 != "nan") {
            sum += $5;
            sum_sq += $5 * $5;
            count++;
            if ($5 > max) max = $5;
            if ($5 < min) min = $5;
        }
    } END {
        mean = sum / count;
        printf "  Std iHS range: %.2f to %.2f\n", min, max;
        printf "  Mean: %.2f\n", mean;
        printf "  Valid scores: %d / %d\n", count, NR;
    }'
    echo ""
else
    echo "ERROR: Output file not created or is empty"
    echo "=========================================="
    exit 1
fi

echo "Next steps:"
echo "  1. Compare EUR-specific iHS with ALL-population iHS"
echo "  2. Validate EUR iHS signals with XP-EHH (EUR vs EAS/AFR)"
echo "  3. Identify true EUR-specific selection signals"
echo ""
