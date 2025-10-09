#!/bin/bash
#SBATCH --job-name=ihs-chr22
#SBATCH --account=mit-ihh-pib
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/%x_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/%x_%j.err

set -euo pipefail
set -x

REPO=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib
HAP=/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin/chr22.hapbin.hap
MAP=/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin/chr22.hapbin.map
OUTDIR=${REPO}/results/ihs
OUTPREFIX=${OUTDIR}/ALL.chr22

# Absolute path to ihsbin (from your earlier command)
IHSBIN=/home/vanbruggenmit/.pixi-detached-envs/mit-ihh-pib-16785450186575381143/envs/default/bin/ihsbin

mkdir -p "$OUTDIR"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Fail early if inputs missing
test -s "$HAP" || { echo "Missing HAP: $HAP"; exit 2; }
test -s "$MAP" || { echo "Missing MAP: $MAP"; exit 2; }

hostname
echo "Using IHSBIN: $IHSBIN with $OMP_NUM_THREADS threads"

# Run
srun "$IHSBIN" \
  --hap "$HAP" \
  --map "$MAP" \
  --out "$OUTPREFIX" \
  --minmaf 0.05 \
  --cutoff 0.1

# Rename to explicit extension for clarity
if [ -s "${OUTPREFIX}" ]; then
  mv -f "${OUTPREFIX}" "${OUTPREFIX}.ihs.tsv"
  echo "Output written to: ${OUTPREFIX}.ihs.tsv"
else
  echo "ERROR: ihsbin produced no output at ${OUTPREFIX}"
  exit 3
fi
