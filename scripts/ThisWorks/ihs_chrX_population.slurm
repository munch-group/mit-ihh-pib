#!/usr/bin/env bash
#SBATCH -J ihs-chrX-POP
#SBATCH -A mit-ihh-pib
#SBATCH -p normal
#SBATCH -c 8
#SBATCH -t 36:00:00
#SBATCH -o /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_chrX_POP_%j.out
#SBATCH -e /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_chrX_POP_%j.err

# Usage: sbatch --export=POP=YRI scripts/ThisWorks/ihs_chrX_population.slurm
# Or: POP=YRI sbatch scripts/ThisWorks/ihs_chrX_population.slurm

set -euo pipefail

# Check if POP is provided
if [[ -z "${POP:-}" ]]; then
    echo "ERROR: POP variable not set"
    echo "Usage: sbatch --export=POP=YRI scripts/ThisWorks/ihs_chrX_population.slurm"
    exit 1
fi

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
DATA_DIR="/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_${POP}"

cd "$REPO_ROOT"

CHR="X"

OUT_DIR="$REPO_ROOT/results/ihs_${POP}"
mkdir -p "$OUT_DIR" "$REPO_ROOT/logs"

HAP="$DATA_DIR/chr${CHR}_${POP}.hapbin.hap"
MAP="$DATA_DIR/chr${CHR}_${POP}.hapbin.map"
OUT="$OUT_DIR/${POP}.chr${CHR}.ihs.tsv"

echo "=========================================="
echo "iHS calculation for chromosome ${CHR}"
echo "Population: ${POP}"
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
    echo "Output already exists: $OUT"
    echo "  If you want to rerun, delete this file first."
    exit 0
fi

# Validate inputs (by claude)
if [[ ! -s "$HAP" ]]; then
    echo "ERROR: HAP file missing or empty: $HAP"
    echo "Did you run: python scripts/ThisWorks/subset_population_hapbin.py ${POP}"
    exit 1
fi

if [[ ! -s "$MAP" ]]; then
    echo "ERROR: MAP file missing or empty: $MAP"
    exit 1
fi

# Check HAP file format (claude)
echo "Verifying HAP file format..."
HAP_LINES=$(wc -l < "$HAP")
HAP_SIZE=$(du -h "$HAP" | cut -f1)
FIELDS_LINE1=$(head -1 "$HAP" | awk '{print NF}')
echo "  HAP file: $HAP_SIZE, $HAP_LINES variants"
echo "  First line: $FIELDS_LINE1 haplotypes (${POP} samples only)"
echo ""

# Check MAP file (claude)
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

echo "Running iHS for chr${CHR} (${POP} samples only)..."
echo "This may take several hours..."
echo ""

pixi run ihsbin --hap "$HAP" --map "$MAP" --out "$OUT"

echo ""
echo "=========================================="
if [[ -s "$OUT" ]]; then
    OUT_LINES=$(wc -l < "$OUT")
    OUT_SIZE=$(du -h "$OUT" | cut -f1)
    echo "iHS completed for chr${CHR} (${POP})"
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
    echo "First statistics:"
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

echo "Analysis complete for ${POP}"
echo ""
