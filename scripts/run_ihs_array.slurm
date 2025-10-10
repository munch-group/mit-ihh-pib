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

# --- absolute locations (as you specified) ---
REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
DATA_HAPBIN_DIR="/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin"

cd "$REPO_ROOT"

CHROMS=(21 22)
CHR="${CHROMS[$SLURM_ARRAY_TASK_ID]}"

THREADS="${THREADS:-8}"

OUT_DIR="$REPO_ROOT/results/ihs"
LOG_DIR="$REPO_ROOT/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

# Inputs (absolute), Outputs (in repo)
HAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.hap"
MAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.map"
OUT="$OUT_DIR/ALL.chr${CHR}.ihs.tsv"

echo "Repo root: $REPO_ROOT"
echo "Inputs:"
echo "  HAP: $HAP"
echo "  MAP: $MAP"
echo "Output:"
echo "  OUT: $OUT"

# Skip if already done
if [[ -s "$OUT" ]]; then
  echo "Already exists: $OUT"
  exit 0
fi

# Validate inputs exist
if [[ ! -s "$HAP" || ! -s "$MAP" ]]; then
  echo "❌ Missing input for chr${CHR}"
  exit 1
fi

echo "Running iHS for chr${CHR} with $THREADS threads..."
pixi run ihsbin "$HAP" "$MAP" "$OUT" "$THREADS"
echo "✅ Done chr${CHR} → $OUT"


