#!/bin/bash
#SBATCH --job-name=test-ihs
#SBATCH --account=mit-ihh-pib
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/%x_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/%x_%j.err
set -euo pipefail
set -x

IHSBIN="/home/vanbruggenmit/.pixi-detached-envs/mit-ihh-pib-16785450186575381143/envs/default/bin/ihsbin"
hostname
srun "$IHSBIN" -h
echo "OK"
