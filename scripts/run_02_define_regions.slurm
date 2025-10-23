#!/usr/bin/env bash
#SBATCH -J ihs_regions
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -o /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_regions_%j.out
#SBATCH -e /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_regions_%j.err

echo "=========================================="
echo "Running iHS Region Definition"
echo "=========================================="
echo "Start time: $(date)"
echo

cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

pixi run python scripts/02_define_selection_regions.py

echo
echo "=========================================="
echo "Region Definition Complete"
echo "=========================================="
echo "End time: $(date)"
