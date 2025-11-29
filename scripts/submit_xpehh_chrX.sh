#!/bin/bash
#SBATCH --job-name=xpehh_chrX
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/xpehh_chrX_%A_%a.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/xpehh_chrX_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-2
#SBATCH --partition=normal

# XP-EHH analysis for chromosome X
# Array job: each task runs one pairwise comparison

set -euo pipefail

# Paths
XPEHH_BIN="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/hapbin_source/build_custom/xpehhbin"
HAP_DIR="/faststorage/project/mit-ihh-pib/data/grch38/hapbin/xpehh"
OUT_DIR="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/results/xp_ehh"

# Create output directory
mkdir -p "${OUT_DIR}"

# Parameters
CUTOFF=0.05       # EHH cutoff value
MINMAF=0.00       # Minimum allele frequency (0 = no filter)
SCALE=20000       # Gap scale parameter in bp
MAX_EXTEND=0      # Maximum distance to traverse (0 = disabled)
BINS=50           # Number of frequency bins for normalization

# Define comparisons
declare -a COMPARISONS=(
    "EUR:AFR"
    "EUR:EAS"
    "AFR:EAS"
)

# Get the comparison for this array task
COMPARISON="${COMPARISONS[$SLURM_ARRAY_TASK_ID]}"
POPA="${COMPARISON%:*}"
POPB="${COMPARISON#*:}"

# Input files
HAP_A="${HAP_DIR}/chrX.${POPA}.hap"
HAP_B="${HAP_DIR}/chrX.${POPB}.hap"
MAP="${HAP_DIR}/chrX.${POPA}.map"
OUT_PREFIX="${OUT_DIR}/chrX_${POPA}_vs_${POPB}"

echo "========================================"
echo "XP-EHH Analysis: ${POPA} vs ${POPB}"
echo "========================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "Date: $(date)"
echo ""
echo "Input files:"
echo "  PopA: ${HAP_A}"
echo "  PopB: ${HAP_B}"
echo "  Map:  ${MAP}"
echo ""
echo "Output: ${OUT_PREFIX}"
echo ""
echo "Parameters:"
echo "  Cutoff: ${CUTOFF}"
echo "  MinMAF: ${MINMAF}"
echo "  Scale: ${SCALE}"
echo "  Max extend: ${MAX_EXTEND}"
echo "  Bins: ${BINS}"
echo "  Threads: ${SLURM_CPUS_PER_TASK}"
echo ""

# Check input files exist
if [[ ! -f "${HAP_A}" ]]; then
    echo "ERROR: ${HAP_A} not found!"
    exit 1
fi
if [[ ! -f "${HAP_B}" ]]; then
    echo "ERROR: ${HAP_B} not found!"
    exit 1
fi
if [[ ! -f "${MAP}" ]]; then
    echo "ERROR: ${MAP} not found!"
    exit 1
fi

# Check if xpehhbin exists
if [[ ! -f "${XPEHH_BIN}" ]]; then
    echo "ERROR: xpehhbin not found at ${XPEHH_BIN}"
    echo "Please build hapbin first:"
    echo "  cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/hapbin_source/build_custom"
    echo "  cmake ../src"
    echo "  make"
    exit 1
fi

# Set OpenMP threads (if hapbin was compiled with OpenMP support)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Run XP-EHH
echo "Starting XP-EHH calculation..."
START_TIME=$(date +%s)

"${XPEHH_BIN}" \
    --hapA "${HAP_A}" \
    --hapB "${HAP_B}" \
    --map "${MAP}" \
    --out "${OUT_PREFIX}" \
    --cutoff ${CUTOFF} \
    --minmaf ${MINMAF} \
    --scale ${SCALE} \
    --max-extend ${MAX_EXTEND} \
    --bin ${BINS}

EXIT_CODE=$?
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo ""
echo "========================================"
if [[ ${EXIT_CODE} -eq 0 ]]; then
    echo "SUCCESS: XP-EHH completed"
else
    echo "ERROR: XP-EHH failed with exit code ${EXIT_CODE}"
fi
echo "========================================"
echo "Duration: ${DURATION} seconds ($(echo "scale=2; ${DURATION}/3600" | bc) hours)"
echo "Output file: ${OUT_PREFIX}"
echo ""

# Show summary statistics if successful
if [[ ${EXIT_CODE} -eq 0 ]] && [[ -f "${OUT_PREFIX}" ]]; then
    echo "Output file statistics:"
    echo "  Total lines: $(wc -l < "${OUT_PREFIX}")"
    NVARIANTS=$(tail -n +2 "${OUT_PREFIX}" | wc -l)
    echo "  Variants analyzed: ${NVARIANTS}"
    echo ""

    # Show first few lines
    echo "First 5 lines of output:"
    head -6 "${OUT_PREFIX}"
    echo ""
fi

echo "Job finished at: $(date)"
exit ${EXIT_CODE}
