#!/usr/bin/env bash
#SBATCH -J ihs_candidates
#SBATCH -A mit-ihh-pib
#SBATCH -p short
#SBATCH -c 1
#SBATCH -t 01:00:00
#SBATCH --mem=16G
#SBATCH -o /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_candidates_%j.out
#SBATCH -e /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_candidates_%j.err

set -euo pipefail

cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

echo "=========================================="
echo "Running iHS Candidate Identification"
echo "=========================================="
echo "Start time: $(date)"
echo ""

# Run the Python analysis script using pixi environment
pixi run python scripts/ThisWorks/01_identify_ihs_candidates.py

echo ""
echo "=========================================="
echo "Analysis Complete!"
echo "End time: $(date)"
echo "=========================================="
