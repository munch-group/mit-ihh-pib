#!/usr/bin/env bash
#SBATCH -J gene_annotation
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -o /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/gene_annotation_%j.out
#SBATCH -e /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/gene_annotation_%j.err

echo "=========================================="
echo "Running Gene Annotation (Phase 3)"
echo "=========================================="
echo "Start time: $(date)"
echo

cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

pixi run python scripts/03_annotate_regions_with_genes.py

echo
echo "=========================================="
echo "Gene Annotation Complete"
echo "=========================================="
echo "End time: $(date)"
