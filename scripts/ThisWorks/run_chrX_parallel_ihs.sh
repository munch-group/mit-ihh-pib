#!/usr/bin/env bash
# Master script to run parallel iHS computation on chromosome X
#
# This splits chromosome X into 5 regions and runs iHS in parallel,
# potentially reducing runtime from 12+ hours to ~3-4 hours
#
# Usage:
#   bash run_chrX_parallel_ihs.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"

cd "$REPO_ROOT"

echo "=========================================="
echo "Parallel iHS Computation for Chromosome X"
echo "=========================================="
echo ""
echo "This will:"
echo "  1. Split chromosome X into 5 genomic regions"
echo "  2. Submit 5 parallel SLURM jobs (6 hours each)"
echo "  3. Merge results when all jobs complete"
echo ""
echo "Expected speedup: ~3-4x faster than sequential"
echo "  (12 hours → ~3-4 hours)"
echo ""

# Step 1: Split the chromosome
echo "Step 1: Splitting chromosome X into regions..."
echo "----------------------------------------"
bash "$SCRIPT_DIR/split_chrX_regions.sh"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to split chromosome X"
    exit 1
fi

echo ""
echo "✓ Chromosome split complete"
echo ""

# Step 2: Submit job array
echo "Step 2: Submitting parallel iHS jobs..."
echo "----------------------------------------"

JOB_ID=$(sbatch --parsable "$SCRIPT_DIR/ihs_chrX_parallel.slurm")

if [[ -z "$JOB_ID" ]]; then
    echo "ERROR: Failed to submit jobs"
    exit 1
fi

echo "✓ Submitted job array: $JOB_ID"
echo ""
echo "Job details:"
echo "  - 5 parallel jobs (regions 1-5)"
echo "  - 8 cores each"
echo "  - 6 hours time limit per job"
echo ""

# Step 3: Submit merge job (depends on array completion)
echo "Step 3: Scheduling merge job..."
echo "----------------------------------------"

cat > "$REPO_ROOT/scripts/ThisWorks/merge_chrX_ihs.slurm" <<'EOF'
#!/usr/bin/env bash
#SBATCH -J merge-chrX-ihs
#SBATCH -A mit-ihh-pib
#SBATCH -p normal
#SBATCH -c 1
#SBATCH -t 1:00:00
#SBATCH -o /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/merge_chrX_ihs_%j.out
#SBATCH -e /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/merge_chrX_ihs_%j.err

bash /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks/merge_chrX_ihs_results.sh
EOF

MERGE_JOB_ID=$(sbatch --parsable --dependency=afterok:${JOB_ID} "$REPO_ROOT/scripts/ThisWorks/merge_chrX_ihs.slurm")

echo "✓ Scheduled merge job: $MERGE_JOB_ID"
echo "  (will run after all region jobs complete)"
echo ""

echo "=========================================="
echo "✓ All jobs submitted!"
echo "=========================================="
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  tail -f $REPO_ROOT/logs/ihs_chrX_region*_${JOB_ID}*.out"
echo ""
echo "Job IDs:"
echo "  Region jobs: ${JOB_ID}_[1-5]"
echo "  Merge job: $MERGE_JOB_ID"
echo ""
echo "Final output will be at:"
echo "  $REPO_ROOT/results/ihs_fixed_par1/ALL.chrX.ihs.tsv"
echo ""
