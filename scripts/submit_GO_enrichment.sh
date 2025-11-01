#!/bin/bash
#SBATCH --job-name=GO_enrichment
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/GO_enrichment_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/GO_enrichment_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal

# GO Enrichment Analysis for X Chromosome Selection Study
# Performs three analyses:
# A. High-confidence genes vs all X genes
# B. High-confidence genes vs X genes under selection
# C. X genes under selection vs autosomal genes under selection

echo "=================================================="
echo "GO Enrichment Analysis"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "=================================================="
echo ""

# Navigate to project directory
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

# Activate pixi environment and run script
/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/.pixi/envs/default/bin/python \
    scripts/run_GO_enrichment_comprehensive.py

exit_code=$?

echo ""
echo "=================================================="
echo "Job completed: $(date)"
echo "Exit code: $exit_code"
echo "=================================================="

exit $exit_code
