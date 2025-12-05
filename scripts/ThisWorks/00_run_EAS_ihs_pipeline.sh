#!/usr/bin/env bash
"""
Master pipeline to run EAS-specific iHS analysis

This pipeline will:
1. Subset chromosome X data to EAS samples only
2. Run iHS on EAS samples
3. Compare EAS iHS with XP-EHH for validation
"""

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
cd "$REPO_ROOT"

echo "=========================================="
echo "EAS-Specific iHS Analysis Pipeline"
echo "=========================================="
echo ""
echo "This pipeline will:"
echo "  1. Subset chrX data to EAS samples (504 samples)"
echo "  2. Run iHS on EAS-only data"
echo "  3. Enable EAS-specific selection signal detection"
echo ""
echo "After EAS iHS completes:"
echo "  → Compare EAS iHS with XP-EHH (EAS vs EUR/AFR)"
echo "  → Identify true EAS-specific selection"
echo "  → Look for known signals (e.g., EDAR pathway)"
echo ""
echo "=========================================="
echo ""

# Step 1: Subset to EAS samples
echo "STEP 1: Subsetting chromosome X to EAS samples..."
echo "----------------------------------------"

if [[ -f "/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_EAS/chrX_EAS.hapbin.hap" ]]; then
    echo "✓ EAS subset already exists"
    echo "  Location: /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_EAS/"
else
    echo "Running subsetting script..."
    python scripts/ThisWorks/01_subset_EAS_samples_hapbin.py

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Subsetting failed"
        exit 1
    fi
fi

echo ""
echo "STEP 2: Submitting iHS job for EAS samples..."
echo "----------------------------------------"

# Check if output already exists
if [[ -f "$REPO_ROOT/results/ihs_EAS/EAS.chrX.ihs.tsv" ]]; then
    echo "✓ EAS iHS output already exists"
    echo "  Location: $REPO_ROOT/results/ihs_EAS/EAS.chrX.ihs.tsv"
    echo ""
    echo "To rerun, delete the output file first."
else
    echo "Submitting SLURM job..."
    JOB_ID=$(sbatch --parsable scripts/ThisWorks/02_ihs_chrX_EAS.slurm)

    echo "✓ Job submitted: $JOB_ID"
    echo ""
    echo "Monitor progress:"
    echo "  squeue -u $USER"
    echo "  tail -f logs/ihs_chrX_EAS_${JOB_ID}.out"
    echo ""
    echo "When complete, run comparison analysis:"
    echo "  python scripts/analyze_EAS_selection.py"
fi

echo ""
echo "=========================================="
echo "Pipeline submitted!"
echo "=========================================="
