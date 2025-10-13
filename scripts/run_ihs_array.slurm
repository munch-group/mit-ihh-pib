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
THREADS="${THREADS:-8}"

OUT_DIR="$REPO_ROOT/results/ihs"
mkdir -p "$OUT_DIR" "$REPO_ROOT/logs"

HAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.hap"
MAP="$DATA_HAPBIN_DIR/chr${CHR}.hapbin.map"
OUT="$OUT_DIR/ALL.chr${CHR}.ihs.tsv"

echo "Repo root: $REPO_ROOT"
echo "Inputs:"
echo "  HAP: $HAP"
echo "  MAP: $MAP"
echo "Output:"
echo "  OUT: $OUT"

[[ -s "$OUT" ]] && { echo "Already exists: $OUT"; exit 0; }
[[ -s "$HAP" && -s "$MAP" ]] || { echo "❌ Missing input for chr${CHR}"; exit 1; }

echo "Running iHS for chr${CHR} with $THREADS threads..."

# Try positional syntax first; if the parser rejects it, retry with long flags.
set +e
pixi run ihsbin "$HAP" "$MAP" "$OUT" "$THREADS"
status=$?
set -e

if [[ $status -ne 0 ]]; then
  echo "Positional args failed (exit $status). Retrying with flags..."
  pixi run ihsbin --hap "$HAP" --map "$MAP" --out "$OUT" --threads "$THREADS"
fi

echo "✅ Done chr${CHR} → $OUT"



