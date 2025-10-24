#!/bin/bash
#SBATCH --job-name=enrichment_analysis
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/enrichment_analysis_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/enrichment_analysis_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=normal

# Phase 4: Functional Enrichment Analysis
# Tests hypothesis: X chromosome genes under selection are enriched for
# reproductive and immune functions

echo "========================================="
echo "Phase 4: Functional Enrichment Analysis"
echo "========================================="
echo "Start time: $(date)"
echo "Running on host: $(hostname)"
echo ""

# Activate conda environment (if you have one)
# Adjust this to your Python environment setup
# source activate your_env_name

# Install gprofiler-official if not already installed
echo "Checking for required packages..."
python3 -c "import gprofiler" 2>/dev/null || {
    echo "Installing gprofiler-official..."
    pip install --user gprofiler-official
}

echo ""
echo "Running enrichment analysis..."
echo ""

# Run the enrichment analysis
python3 /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/04_functional_enrichment_analysis.py

exit_code=$?

echo ""
echo "========================================="
echo "Job completed with exit code: $exit_code"
echo "End time: $(date)"
echo "========================================="

exit $exit_code
