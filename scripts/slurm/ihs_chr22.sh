#!/bin/bash
#SBATCH --job-name=ihs-chr22
#SBATCH --account=mit-ihh-pib
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/%x_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/%x_%j.err
set -euo pipefail

# 1) Go to the folder that has pyproject.toml
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# 2) Run hapbin iHS
pixi exec ihsbin \
  --hap /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin/chr22.hapbin.hap \
  --map /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin/chr22.hapbin.map \
  --out /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/ihs/ALL.chr22 \
  --minmaf 0.05 \
  --cutoff 0.1
