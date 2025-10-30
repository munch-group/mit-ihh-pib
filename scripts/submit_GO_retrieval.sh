#!/bin/bash
#SBATCH --job-name=GO_retrieval_X
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/GO_retrieval_X_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/GO_retrieval_X_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

# GO Term Retrieval for X Chromosome Genes
# This job retrieves GO annotations for ~664 X chromosome genes
# Estimated runtime: 6-8 hours

echo "=================================================="
echo "GO Term Retrieval for X Chromosome Genes"
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "=================================================="
echo ""

# Navigate to project directory
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

# Activate pixi environment and run script
/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/.pixi/envs/default/bin/python \
    scripts/retrieve_GO_terms_X_genes.py

exit_code=$?

echo ""
echo "=================================================="
echo "Job completed: $(date)"
echo "Exit code: $exit_code"
echo "=================================================="

exit $exit_code
