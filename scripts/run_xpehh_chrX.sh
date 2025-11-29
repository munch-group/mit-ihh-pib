#!/bin/bash
#
# Run XP-EHH analysis for chromosome X
# Compares all pairwise combinations of populations
#

set -euo pipefail

# Paths
XPEHH_BIN="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/hapbin_source/build_custom/xpehhbin"
HAP_DIR="/faststorage/project/mit-ihh-pib/data/grch38/hapbin/xpehh"
OUT_DIR="/faststorage/project/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/xpehh/chrX"

# Create output directory
mkdir -p "${OUT_DIR}"

# Parameters (adjust as needed)
CUTOFF=0.05       # EHH cutoff value
MINMAF=0.00       # Minimum allele frequency
SCALE=20000       # Gap scale parameter in bp
MAX_EXTEND=0      # Maximum distance to traverse (0 = disabled)
BINS=50           # Number of frequency bins for normalization

# Population codes
POPS=("EUR" "AFR" "EAS")

echo "========================================"
echo "XP-EHH Analysis for Chromosome X"
echo "========================================"
echo "Input directory: ${HAP_DIR}"
echo "Output directory: ${OUT_DIR}"
echo ""
echo "Parameters:"
echo "  Cutoff: ${CUTOFF}"
echo "  MinMAF: ${MINMAF}"
echo "  Scale: ${SCALE}"
echo "  Max extend: ${MAX_EXTEND}"
echo "  Bins: ${BINS}"
echo ""

# Run all pairwise comparisons
for ((i=0; i<${#POPS[@]}; i++)); do
    for ((j=i+1; j<${#POPS[@]}; j++)); do
        POPA="${POPS[i]}"
        POPB="${POPS[j]}"

        HAP_A="${HAP_DIR}/chrX.${POPA}.hap"
        HAP_B="${HAP_DIR}/chrX.${POPB}.hap"
        MAP="${HAP_DIR}/chrX.${POPA}.map"
        OUT_PREFIX="${OUT_DIR}/chrX_${POPA}_vs_${POPB}"

        echo "========================================"
        echo "Running: ${POPA} vs ${POPB}"
        echo "========================================"
        echo "PopA: ${HAP_A}"
        echo "PopB: ${HAP_B}"
        echo "Map:  ${MAP}"
        echo "Output: ${OUT_PREFIX}"
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

        END_TIME=$(date +%s)
        DURATION=$((END_TIME - START_TIME))

        echo ""
        echo "Completed in ${DURATION} seconds"
        echo "Output file: ${OUT_PREFIX}"
        echo ""

        # Show summary statistics
        if [[ -f "${OUT_PREFIX}" ]]; then
            NVARIANTS=$(tail -n +2 "${OUT_PREFIX}" | wc -l)
            echo "Number of variants analyzed: ${NVARIANTS}"
            echo ""
        fi
    done
done

echo "========================================"
echo "All XP-EHH analyses completed!"
echo "========================================"
echo ""
echo "Output files:"
ls -lh "${OUT_DIR}"/chrX_*

echo ""
echo "Done!"
