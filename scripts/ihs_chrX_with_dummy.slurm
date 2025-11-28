#!/usr/bin/env bash
#SBATCH -J ihs-chrX-dummy
#SBATCH -A mit-ihh-pib
#SBATCH -p normal
#SBATCH -c 8
#SBATCH -t 36:00:00
#SBATCH -o /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_chrX_dummy_%j.out
#SBATCH -e /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_chrX_dummy_%j.err

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
DATA_HAPBIN_DIR="/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_with_dummy"

# Use modified ihsbin that supports '-' dummy characters
IHSBIN="$REPO_ROOT/hapbin_source/build_custom/ihsbin"

cd "$REPO_ROOT"

CHR="X"

OUT_DIR="$REPO_ROOT/results/ihs_with_dummy"
mkdir -p "$OUT_DIR" "$REPO_ROOT/logs"

HAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.hap"
MAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.map"
OUT="$OUT_DIR/ALL.chr${CHR}.ihs.tsv"

echo "=========================================="
echo "iHS calculation for chromosome ${CHR}"
echo "Using modified ihsbin with dummy '-' support"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================="
echo "Inputs:"
echo "  HAP: $HAP"
echo "  MAP: $MAP"
echo "  IHSBIN: $IHSBIN"
echo "Output:"
echo "  OUT: $OUT"
echo ""

# Verify modified ihsbin exists
if [[ ! -x "$IHSBIN" ]]; then
    echo "ERROR: Modified ihsbin binary not found or not executable: $IHSBIN"
    exit 1
fi

# Validate inputs exist (skip slow verification on 31GB files)
if [[ ! -s "$HAP" ]]; then
    echo "ERROR: HAP file missing or empty: $HAP"
    exit 1
fi

if [[ ! -s "$MAP" ]]; then
    echo "ERROR: MAP file missing or empty: $MAP"
    exit 1
fi

echo "Input files verified (skipping slow checks on large files)"
echo ""

echo "Running iHS for chr${CHR} with modified ihsbin..."
echo ""

"$IHSBIN" --hap "$HAP" --map "$MAP" --out "$OUT"

echo ""
echo "=========================================="
if [[ -s "$OUT" ]]; then
    OUT_SIZE=$(du -h "$OUT" | cut -f1)
    echo "✓ iHS completed for chr${CHR}"
    echo "=========================================="
    echo ""
    echo "Output file: $OUT"
    echo "  Size: $OUT_SIZE"
    echo ""

    # Show first few results
    echo "First 5 iHS results:"
    head -6 "$OUT"
    echo ""

    # Check if PAR1 is included
    echo "Checking PAR1 coverage (expect variants starting ~288,060):"
    tail -n +2 "$OUT" | head -5 | awk '{print "  Position:", $2}'
    echo ""
else
    echo "ERROR: Output file not created or is empty"
    echo "=========================================="
    exit 1
fi

echo "Expected: iHS output should include PAR1 variants with proper male hemizygosity handling!"
echo ""
