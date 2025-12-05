#!/usr/bin/env bash
"""
Master pipeline to run AFR-specific iHS analysis

This pipeline will:
1. Subset chromosome X data to AFR samples only
2. Run iHS on AFR samples
3. Compare AFR iHS with XP-EHH for validation
"""

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
cd "$REPO_ROOT"

echo "=========================================="
echo "AFR-Specific iHS Analysis Pipeline"
echo "=========================================="
echo ""
echo "This pipeline will:"
echo "  1. Subset chrX data to AFR samples (504 samples)"
echo "  2. Run iHS on AFR-only data"
echo "  3. Enable AFR-specific selection signal detection"
echo ""
echo "After AFR iHS completes:"
echo "  → Compare AFR iHS with XP-EHH (AFR vs EAS/EUR)"
echo "  → Identify true AFR-specific selection"
echo ""
echo "=========================================="
echo ""

# Step 1: Subset to AFR samples
echo "STEP 1: Subsetting chromosome X to AFR samples..."
echo "----------------------------------------"

if [[ -f "/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_AFR/chrX_AFR.hapbin.hap" ]]; then
    echo "✓ AFR subset already exists"
    echo "  Location: /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_AFR/"
else
    echo "Running subsetting script..."
    python scripts/ThisWorks/01_subset_AFR_samples_hapbin.py

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Subsetting failed"
        exit 1
    fi
fi

echo ""
echo "STEP 2: Submitting iHS job for AFR samples..."
echo "----------------------------------------"

# Check if output already exists
if [[ -f "$REPO_ROOT/results/ihs_AFR/AFR.chrX.ihs.tsv" ]]; then
    echo "✓ AFR iHS output already exists"
    echo "  Location: $REPO_ROOT/results/ihs_AFR/AFR.chrX.ihs.tsv"
    echo ""
    echo "To rerun, delete the output file first."
else
    echo "Submitting SLURM job..."
    JOB_ID=$(sbatch --parsable scripts/ThisWorks/02_ihs_chrX_AFR.slurm)

    echo "✓ Job submitted: $JOB_ID"
    echo ""
    echo "Monitor progress:"
    echo "  squeue -u $USER"
    echo "  tail -f logs/ihs_chrX_AFR_${JOB_ID}.out"
    echo ""
    echo "When complete, run comparison analysis:"
    echo "  python scripts/analyze_AFR_selection.py"
fi

echo ""
echo "=========================================="
echo "Pipeline submitted!"
echo "=========================================="
