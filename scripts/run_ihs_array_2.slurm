#!/usr/bin/env bash
#SBATCH -J ihs-21-22
#SBATCH -A mit-ihh-pib
#SBATCH -p normal
#SBATCH -c 8
#SBATCH -t 02:00:00
#SBATCH -o /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs-%A_%a.out
#SBATCH -e /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs-%A_%a.err
#SBATCH --array=0-1

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
DATA_HAPBIN_DIR="/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin"

cd "$REPO_ROOT"

CHROMS=(21 22)
CHR="${CHROMS[$SLURM_ARRAY_TASK_ID]}"

OUT_DIR="$REPO_ROOT/results/ihs"
mkdir -p "$OUT_DIR" "$REPO_ROOT/logs"

HAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.hap"
MAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.map"
OUT="$OUT_DIR/ALL.chr${CHR}.ihs.tsv"

echo "=========================================="
echo "iHS calculation for chromosome ${CHR}"
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
    exit 0
fi

# Validate inputs
if [[ ! -s "$HAP" ]]; then
    echo "ERROR: HAP file missing or empty: $HAP"
    exit 1
fi

if [[ ! -s "$MAP" ]]; then
    echo "ERROR: MAP file missing or empty: $MAP"
    exit 1
fi

echo "Running iHS for chr${CHR}..."

# ihsbin only accepts --hap, --map, and --out flags
# It does NOT accept --threads
pixi run ihsbin --hap "$HAP" --map "$MAP" --out "$OUT"

if [[ -s "$OUT" ]]; then
    echo ""
    echo "Done chr${CHR} → $OUT"
    echo "   Output file size: $(du -h "$OUT" | cut -f1)"
else
    echo "ERROR: Output file not created or is empty"
    exit 1
fi