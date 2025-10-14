#!/usr/bin/env bash
#SBATCH --job-name=hapbin_chrX
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/hapbin_chrX_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/hapbin_chrX_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal

# Set chromosome to X
CHR="X"

DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"
SCRIPT_PATH="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks/convert_to_hapbin.py"
OUT_DIR="$DATA_ROOT/hapbin"

mkdir -p "$OUT_DIR"

echo "=========================================="
echo "Processing chromosome $CHR (Job ID: $SLURM_JOB_ID)"
echo "=========================================="

# Find HAP file
HAP=""
LEG=""
if [[ -f "$DATA_ROOT/work/chr${CHR}/chr${CHR}.impute.hap.gz" ]]; then
    HAP="$DATA_ROOT/work/chr${CHR}/chr${CHR}.impute.hap.gz"
    LEG="$DATA_ROOT/work/chr${CHR}/chr${CHR}.impute.legend.gz"
    echo "Found HAP/LEGEND in: work/chr${CHR}/"
elif [[ -f "$DATA_ROOT/raw/chr${CHR}/chr${CHR}.impute.hap.gz" ]]; then
    HAP="$DATA_ROOT/raw/chr${CHR}/chr${CHR}.impute.hap.gz"
    LEG="$DATA_ROOT/raw/chr${CHR}/chr${CHR}.impute.legend.gz"
    echo "Found HAP/LEGEND in: raw/chr${CHR}/"
elif [[ -f "$DATA_ROOT/impute/chr${CHR}.hap.gz" ]]; then
    HAP="$DATA_ROOT/impute/chr${CHR}.hap.gz"
    LEG="$DATA_ROOT/impute/chr${CHR}.legend.gz"
    echo "Found HAP/LEGEND in: impute/"
fi

if [[ -z "$HAP" ]]; then
    echo "  ✗ HAP file not found for chr$CHR"
    exit 1
fi

ANC="$DATA_ROOT/raw/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${CHR}.fa"
# Use the adjusted recombination map for chrX (with 2/3 multiplier)
MAP="$DATA_ROOT/work/chr${CHR}/genetic_map_chrX_adjusted.txt"

# Fallback to standard map location if adjusted one doesn't exist
if [[ ! -f "$MAP" ]]; then
    echo "  ⚠ Adjusted map not found, checking standard location..."
    MAP="$DATA_ROOT/maps/chr${CHR}.decode.sexavg.cm.tsv"
fi

if [[ ! -f "$ANC" ]]; then
    echo "  ✗ Ancestral FASTA not found: $ANC"
    exit 1
fi

if [[ ! -f "$MAP" ]]; then
    echo "  ✗ Recombination map not found: $MAP"
    exit 1
fi

echo "Using recombination map: $MAP"

# Auto-detect chromosome name
if [[ -f "$ANC.fai" ]]; then
    FASTA_CHR=$(head -1 "$ANC.fai" | awk '{print $1}')
else
    FASTA_CHR=$(grep "^>" "$ANC" | head -1 | sed 's/^>//' | awk '{print $1}')
fi

echo "Using FASTA chromosome name: $FASTA_CHR"
echo "Using Python script: $SCRIPT_PATH"
echo ""

# Verify script exists
if [[ ! -f "$SCRIPT_PATH" ]]; then
    echo "ERROR: Python script not found at $SCRIPT_PATH"
    exit 1
fi

# Run conversion
python3 "$SCRIPT_PATH" \
    --hap "$HAP" \
    --legend "$LEG" \
    --anc-fasta "$ANC" \
    --chr "$FASTA_CHR" \
    --recomb-tsv "$MAP" \
    --no-header \
    --pos-col 1 \
    --cm-col 2 \
    --out-prefix "$OUT_DIR/chr${CHR}.hapbin" \
    --verbose

if [ $? -eq 0 ]; then
    echo ""
    echo "  ✓ Successfully converted chr$CHR"
    if [[ -f "$OUT_DIR/chr${CHR}.hapbin.hap" ]]; then
        HAP_SIZE=$(du -h "$OUT_DIR/chr${CHR}.hapbin.hap" | cut -f1)
        HAP_LINES=$(wc -l < "$OUT_DIR/chr${CHR}.hapbin.hap")
        echo "  Output: $HAP_LINES SNPs ($HAP_SIZE)"
    fi
else
    echo ""
    echo "  ✗ Failed to convert chr$CHR"
    exit 1
fi

echo ""
echo "=========================================="
echo "Chromosome X processing complete"
echo "=========================================="