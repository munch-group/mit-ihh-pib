#!/usr/bin/env bash
"""
Master pipeline to run population-specific iHS analysis for chromosomes 1-22

This pipeline will:
1. Subset chromosomes 1-22 data to AFR, EUR, and EAS samples
2. Run iHS on all chromosomes for each population using SLURM job arrays
3. Monitor progress and provide status updates

Populations: AFR (African), EUR (European), EAS (East Asian)
Chromosomes: chr1-chr22 (autosomes only)
"""

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
cd "$REPO_ROOT"

echo "=========================================="
echo "Population-Specific iHS Analysis Pipeline"
echo "=========================================="
echo ""
echo "This pipeline will:"
echo "  1. Subset chr1-22 data to AFR, EUR, and EAS samples"
echo "  2. Run iHS on 22 chromosomes × 3 populations = 66 jobs"
echo "  3. Store results in population-specific directories"
echo ""
echo "Populations:"
echo "  - AFR: African ancestry"
echo "  - EUR: European ancestry"
echo "  - EAS: East Asian ancestry"
echo ""
echo "Output directories:"
echo "  - results/ihs_AFR/"
echo "  - results/ihs_EUR/"
echo "  - results/ihs_EAS/"
echo ""
echo "Log directory:"
echo "  - logs/"
echo ""
echo "=========================================="
echo ""

# Create necessary directories
mkdir -p "$REPO_ROOT/logs"
mkdir -p "$REPO_ROOT/results/ihs_AFR"
mkdir -p "$REPO_ROOT/results/ihs_EUR"
mkdir -p "$REPO_ROOT/results/ihs_EAS"

# Step 1: Subset to population-specific samples
echo "STEP 1: Subsetting chromosomes 1-22 to AFR, EUR, and EAS samples..."
echo "----------------------------------------"
echo ""

# Check if subsetting is already complete by checking a few key files
AFR_CHR1="/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_AFR/chr1_AFR.hapbin.hap"
EUR_CHR1="/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_EUR/chr1_EUR.hapbin.hap"
EAS_CHR1="/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_EAS/chr1_EAS.hapbin.hap"

if [[ -f "$AFR_CHR1" ]] && [[ -f "$EUR_CHR1" ]] && [[ -f "$EAS_CHR1" ]]; then
    echo "✓ Population subsets appear to exist"
    echo "  Checking completeness..."
    echo ""

    # Count existing files for each population
    AFR_COUNT=$(ls /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_AFR/chr*_AFR.hapbin.hap 2>/dev/null | wc -l)
    EUR_COUNT=$(ls /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_EUR/chr*_EUR.hapbin.hap 2>/dev/null | wc -l)
    EAS_COUNT=$(ls /home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_EAS/chr*_EAS.hapbin.hap 2>/dev/null | wc -l)

    echo "  AFR: $AFR_COUNT/22 chromosomes"
    echo "  EUR: $EUR_COUNT/22 chromosomes"
    echo "  EAS: $EAS_COUNT/22 chromosomes"
    echo ""

    if [[ $AFR_COUNT -eq 22 ]] && [[ $EUR_COUNT -eq 22 ]] && [[ $EAS_COUNT -eq 22 ]]; then
        echo "✓ All population subsets complete (22 chromosomes × 3 populations)"
        echo ""
    else
        echo "⚠ Subsetting incomplete. Running subsetting script..."
        echo ""
        python scripts/ThisWorks/01_subset_populations_chr1-22.py

        if [[ $? -ne 0 ]]; then
            echo "ERROR: Subsetting failed"
            exit 1
        fi
    fi
else
    echo "Running subsetting script for the first time..."
    echo "This will create population-specific .hap and .map files for chr1-22"
    echo ""
    python scripts/ThisWorks/01_subset_populations_chr1-22.py

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Subsetting failed"
        exit 1
    fi
fi

echo ""
echo "STEP 2: Submitting iHS job array..."
echo "----------------------------------------"
echo ""

# Check how many jobs are already complete
COMPLETED_JOBS=0
TOTAL_JOBS=66

for POP in AFR EUR EAS; do
    for CHR in {1..22}; do
        OUT_FILE="$REPO_ROOT/results/ihs_${POP}/${POP}.chr${CHR}.ihs.tsv"
        if [[ -s "$OUT_FILE" ]]; then
            ((COMPLETED_JOBS++))
        fi
    done
done

echo "Current progress: $COMPLETED_JOBS/$TOTAL_JOBS jobs completed"
echo ""

if [[ $COMPLETED_JOBS -eq $TOTAL_JOBS ]]; then
    echo "✓ All iHS jobs already completed!"
    echo ""
    echo "Output locations:"
    echo "  - $REPO_ROOT/results/ihs_AFR/"
    echo "  - $REPO_ROOT/results/ihs_EUR/"
    echo "  - $REPO_ROOT/results/ihs_EAS/"
    echo ""
    echo "To rerun specific jobs, delete the corresponding output files."
else
    echo "Submitting SLURM job array (66 jobs: 22 chromosomes × 3 populations)..."
    echo ""

    JOB_ID=$(sbatch --parsable scripts/ThisWorks/02_ihs_chr1-22_populations.slurm)

    echo "✓ Job array submitted!"
    echo "=========================================="
    echo ""
    echo "Job Array ID: $JOB_ID"
    echo "Total jobs: 66 (22 chromosomes × 3 populations)"
    echo ""
    echo "Monitor progress:"
    echo "  squeue -u $USER"
    echo "  squeue -j $JOB_ID"
    echo ""
    echo "Check individual job logs:"
    echo "  ls -lth logs/ihs_autosomes_${JOB_ID}_*.out | head"
    echo "  tail -f logs/ihs_autosomes_${JOB_ID}_1.out"
    echo ""
    echo "Check for errors:"
    echo "  grep -l ERROR logs/ihs_autosomes_${JOB_ID}_*.err"
    echo ""
    echo "Count completed jobs:"
    echo "  ls results/ihs_AFR/*.ihs.tsv 2>/dev/null | wc -l"
    echo "  ls results/ihs_EUR/*.ihs.tsv 2>/dev/null | wc -l"
    echo "  ls results/ihs_EAS/*.ihs.tsv 2>/dev/null | wc -l"
    echo ""
    echo "Output directories:"
    echo "  AFR: $REPO_ROOT/results/ihs_AFR/"
    echo "  EUR: $REPO_ROOT/results/ihs_EUR/"
    echo "  EAS: $REPO_ROOT/results/ihs_EAS/"
    echo ""
fi

echo "=========================================="
echo "Pipeline setup complete!"
echo "=========================================="
echo ""
echo "Job breakdown:"
echo "  Task 1-22:   AFR chr1-chr22"
echo "  Task 23-44:  EUR chr1-chr22"
echo "  Task 45-66:  EAS chr1-chr22"
echo ""
echo "Expected runtime: ~24-48 hours (depends on cluster load)"
echo ""
