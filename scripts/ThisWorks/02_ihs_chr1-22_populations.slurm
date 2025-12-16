#!/usr/bin/env bash
#SBATCH -J ihs-autosomes
#SBATCH -A mit-ihh-pib
#SBATCH -p normal
#SBATCH -c 8
#SBATCH -t 48:00:00
#SBATCH --array=1-66
#SBATCH -o /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_autosomes_%A_%a.out
#SBATCH -e /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_autosomes_%A_%a.err

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
cd "$REPO_ROOT"

# Create logs directory if it doesn't exist
mkdir -p "$REPO_ROOT/logs"

# Define populations and chromosomes
POPULATIONS=("AFR" "EUR" "EAS")
CHROMOSOMES=($(seq 1 22))

# Calculate which population and chromosome to process
# Array index goes from 1 to 66 (22 chromosomes × 3 populations)
TASK_ID=$((SLURM_ARRAY_TASK_ID - 1))
POP_IDX=$((TASK_ID / 22))
CHR_IDX=$((TASK_ID % 22))

POP="${POPULATIONS[$POP_IDX]}"
CHR="${CHROMOSOMES[$CHR_IDX]}"

echo "=========================================="
echo "iHS calculation"
echo "Job Array ID: $SLURM_ARRAY_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "=========================================="
echo "Population: $POP"
echo "Chromosome: chr${CHR}"
echo "=========================================="
echo ""

# Set paths
DATA_DIR="/home/vanbruggenmit/mit-ihh-pib/data/grch38/hapbin_${POP}"
OUT_DIR="$REPO_ROOT/results/ihs_${POP}"
mkdir -p "$OUT_DIR"

HAP="$DATA_DIR/chr${CHR}_${POP}.hapbin.hap"
MAP="$DATA_DIR/chr${CHR}_${POP}.hapbin.map"
OUT="$OUT_DIR/${POP}.chr${CHR}.ihs.tsv"

echo "Inputs:"
echo "  HAP: $HAP"
echo "  MAP: $MAP"
echo "Output:"
echo "  OUT: $OUT"
echo ""

# Check if output already exists
if [[ -s "$OUT" ]]; then
    echo "✓ Output already exists: $OUT"
    echo "  If you want to rerun, delete this file first."
    exit 0
fi

# Validate inputs
if [[ ! -s "$HAP" ]]; then
    echo "ERROR: HAP file missing or empty: $HAP"
    echo "Did you run 01_subset_populations_chr1-22.py first?"
    exit 1
fi

if [[ ! -s "$MAP" ]]; then
    echo "ERROR: MAP file missing or empty: $MAP"
    exit 1
fi

# Check HAP file format
echo "Verifying HAP file format..."
HAP_LINES=$(wc -l < "$HAP")
HAP_SIZE=$(du -h "$HAP" | cut -f1)
FIELDS_LINE1=$(head -1 "$HAP" | awk '{print NF}')
echo "  HAP file: $HAP_SIZE, $HAP_LINES variants"
echo "  Haplotypes: $FIELDS_LINE1 (${POP} samples only)"
echo ""

# Check MAP file
echo "Verifying MAP file..."
MAP_LINES=$(wc -l < "$MAP")
MAP_FIRST_POS=$(head -1 "$MAP" | awk '{print $4}')
MAP_LAST_POS=$(tail -1 "$MAP" | awk '{print $4}')
echo "  MAP file: $MAP_LINES positions"
echo "  Position range: $MAP_FIRST_POS - $MAP_LAST_POS"
echo ""

if [[ "$HAP_LINES" != "$MAP_LINES" ]]; then
    echo "WARNING: HAP and MAP have different line counts"
    echo "  HAP: $HAP_LINES"
    echo "  MAP: $MAP_LINES"
fi

echo "Running iHS for chr${CHR} (${POP} samples only)..."
echo "This may take several hours..."
echo ""

# Record start time
START_TIME=$(date +%s)

pixi run ihsbin --hap "$HAP" --map "$MAP" --out "$OUT"

# Record end time
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))

echo ""
echo "=========================================="
if [[ -s "$OUT" ]]; then
    OUT_LINES=$(wc -l < "$OUT")
    OUT_SIZE=$(du -h "$OUT" | cut -f1)
    echo "✓ iHS completed for chr${CHR} (${POP})"
    echo "=========================================="
    echo ""
    echo "Output file: $OUT"
    echo "  Size: $OUT_SIZE"
    echo "  Lines: $OUT_LINES"
    echo "  Runtime: ${HOURS}h ${MINUTES}m"
    echo ""

    # Show first few results
    echo "First 5 iHS results:"
    head -6 "$OUT"
    echo ""

    # Statistics
    echo "Quick statistics:"
    tail -n +2 "$OUT" | awk '{
        if ($5 != "NA" && $5 != "nan") {
            sum += $5;
            sum_sq += $5 * $5;
            count++;
            if (NR == 2 || $5 > max) max = $5;
            if (NR == 2 || $5 < min) min = $5;
        }
    } END {
        if (count > 0) {
            mean = sum / count;
            printf "  Std iHS range: %.2f to %.2f\n", min, max;
            printf "  Mean: %.2f\n", mean;
            printf "  Valid scores: %d / %d\n", count, NR;
        } else {
            printf "  No valid iHS scores found\n";
        }
    }'
    echo ""
else
    echo "ERROR: Output file not created or is empty"
    echo "=========================================="
    exit 1
fi

echo "Job completed successfully!"
echo ""
