#!/usr/bin/env bash
#SBATCH -J ihs-21-22
#SBATCH -A mit-ihh-pib
#SBATCH -p normal
#SBATCH -c 8
#SBATCH -t 02:00:00
#SBATCH -o logs/ihs-%A_%a.out
#SBATCH -e logs/ihs-%A_%a.err
#SBATCH --array=0-1

set -euo pipefail
REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"

CHROMS=(21 22)
CHR="${CHROMS[$SLURM_ARRAY_TASK_ID]}"

THREADS="${THREADS:-8}"
HAP_DIR="results/hapbin"
MAP_DIR="results/hapbin"
OUT_DIR="results/ihs"
mkdir -p "$OUT_DIR" logs

HAP="${HAP_DIR}/ALL.chr${CHR}.hap"
MAP="${MAP_DIR}/ALL.chr${CHR}.map"
OUT="${OUT_DIR}/ALL.chr${CHR}.ihs.tsv"

[[ -s "$OUT" ]] && { echo "Already exists: $OUT"; exit 0; }
[[ -s "$HAP" && -s "$MAP" ]] || { echo "Missing input for chr${CHR}"; exit 1; }

echo "Running iHS chr${CHR} with $THREADS threads"
pixi run ihsbin "$HAP" "$MAP" "$OUT" "$THREADS"
echo "Done chr${CHR} → $OUT"

