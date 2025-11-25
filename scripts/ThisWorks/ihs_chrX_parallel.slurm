#!/usr/bin/env bash
#SBATCH -J ihs-chrX-region
#SBATCH -A mit-ihh-pib
#SBATCH -p normal
#SBATCH -c 8
#SBATCH -t 6:00:00
#SBATCH --array=1-5
#SBATCH -o /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_chrX_region%a_%j.out
#SBATCH -e /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_chrX_region%a_%j.err

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
DATA_REGION_DIR="/faststorage/project/mit-ihh-pib/data/grch38/hapbin/chrX_regions"

cd "$REPO_ROOT"

# Get region number from SLURM array task ID
REGION_NUM=$SLURM_ARRAY_TASK_ID

OUT_DIR="$REPO_ROOT/results/ihs_fixed_par1_paraproc/regions"
mkdir -p "$OUT_DIR" "$REPO_ROOT/logs"

HAP="$DATA_REGION_DIR/chrX_region${REGION_NUM}.hap"
MAP="$DATA_REGION_DIR/chrX_region${REGION_NUM}.map"
OUT="$OUT_DIR/chrX_region${REGION_NUM}.ihs.tsv"

echo "=========================================="
echo "iHS calculation for chromosome X - Region ${REGION_NUM}"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "=========================================="
echo "Inputs:"
echo "  HAP: $HAP"
echo "  MAP: $MAP"
echo "Output:"
echo "  OUT: $OUT"
echo ""

# Check HAP file format
if [[ ! -s "$HAP" ]]; then
    echo "ERROR: HAP file missing or empty: $HAP"
    echo "Have you run split_chrX_regions.sh first?"
    exit 1
fi

if [[ ! -s "$MAP" ]]; then
    echo "ERROR: MAP file missing or empty: $MAP"
    echo "Have you run split_chrX_regions.sh first?"
    exit 1
fi

echo "Verifying input files..."
HAP_LINES=$(wc -l < "$HAP")
MAP_LINES=$(wc -l < "$MAP")
HAP_SIZE=$(du -h "$HAP" | cut -f1)
FIELDS_LINE1=$(head -1 "$HAP" | awk '{print NF}')
echo "  HAP file: $HAP_SIZE, $HAP_LINES SNPs, $FIELDS_LINE1 haplotypes"
echo "  MAP file: $MAP_LINES positions"

MAP_FIRST_POS=$(head -1 "$MAP" | awk '{print $4}')
MAP_LAST_POS=$(tail -1 "$MAP" | awk '{print $4}')
echo "  Position range: $MAP_FIRST_POS - $MAP_LAST_POS"
echo ""

if [[ "$HAP_LINES" != "$MAP_LINES" ]]; then
    echo "ERROR: HAP and MAP have different line counts"
    echo "  HAP: $HAP_LINES"
    echo "  MAP: $MAP_LINES"
    exit 1
fi

echo "Running iHS for region ${REGION_NUM}..."
echo ""

pixi run ihsbin --hap "$HAP" --map "$MAP" --out "$OUT"

echo ""
echo "=========================================="
if [[ -s "$OUT" ]]; then
    OUT_LINES=$(wc -l < "$OUT")
    OUT_SIZE=$(du -h "$OUT" | cut -f1)
    echo "✓ iHS completed for region ${REGION_NUM}"
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
else
    echo "ERROR: Output file not created or is empty"
    echo "=========================================="
    exit 1
fi

echo "Region ${REGION_NUM} complete!"
echo ""
