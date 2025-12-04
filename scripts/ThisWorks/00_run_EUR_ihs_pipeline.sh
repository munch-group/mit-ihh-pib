#!/usr/bin/env bash
"""
Master pipeline to run EUR-specific iHS analysis

This pipeline will:
1. Subset chromosome X data to EUR samples only
2. Run iHS on EUR samples
3. Compare EUR iHS with XP-EHH for validation
"""

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
cd "$REPO_ROOT"

echo "=========================================="
echo "EUR-Specific iHS Analysis Pipeline"
echo "=========================================="
echo ""
echo "This pipeline will:"
echo "  1. Subset chrX data to EUR samples (503 samples)"
echo "  2. Run iHS on EUR-only data"
echo "  3. Enable EUR-specific selection signal detection"
echo ""
echo "Current analysis uses ALL populations (1,668 samples)"
echo "  → Cannot distinguish EUR-specific from shared signals"
echo ""
echo "After EUR iHS completes:"
echo "  → Compare EUR iHS with XP-EHH (EUR vs EAS/AFR)"
echo "  → Identify true EUR-specific selection"
echo ""
echo "=========================================="
echo ""

# Step 1: Subset to EUR samples
echo "STEP 1: Subsetting chromosome X to EUR samples..."
echo "----------------------------------------"

if [[ -f "/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_EUR/chrX_EUR.hapbin.hap" ]]; then
    echo "✓ EUR subset already exists"
    echo "  Location: /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_EUR/"
else
    echo "Running subsetting script..."
    python scripts/ThisWorks/01_subset_EUR_samples_hapbin.py

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Subsetting failed"
        exit 1
    fi
fi

echo ""
echo "STEP 2: Submitting iHS job for EUR samples..."
echo "----------------------------------------"

# Check if output already exists
if [[ -f "$REPO_ROOT/results/ihs_EUR/EUR.chrX.ihs.tsv" ]]; then
    echo "✓ EUR iHS output already exists"
    echo "  Location: $REPO_ROOT/results/ihs_EUR/EUR.chrX.ihs.tsv"
    echo ""
    echo "To rerun, delete the output file first."
else
    echo "Submitting SLURM job..."
    JOB_ID=$(sbatch --parsable scripts/ThisWorks/02_ihs_chrX_EUR.slurm)

    echo "✓ Job submitted: $JOB_ID"
    echo ""
    echo "Monitor progress:"
    echo "  squeue -u $USER"
    echo "  tail -f logs/ihs_chrX_EUR_${JOB_ID}.out"
    echo ""
    echo "When complete, run comparison analysis:"
    echo "  python scripts/ThisWorks/03_compare_EUR_vs_ALL_ihs.py"
fi

echo ""
echo "=========================================="
echo "Pipeline submitted!"
echo "=========================================="
