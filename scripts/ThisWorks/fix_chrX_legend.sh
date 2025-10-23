#!/usr/bin/env bash
#SBATCH --job-name=fix_chrX_legend
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/fix_chrX_legend_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/fix_chrX_legend_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal

set -euo pipefail

CHR="X"
DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"
RAW_DIR="$DATA_ROOT/raw/chr${CHR}"
WORK_DIR="$DATA_ROOT/work/chr${CHR}"

VCF_FILE="$RAW_DIR/chr${CHR}.vcf.gz"
LEGEND_OUT="$WORK_DIR/chr${CHR}.impute.legend"

echo "=========================================="
echo "Fixing chrX legend file with proper IDs"
echo "=========================================="
echo ""

# Check if VCF exists
if [[ ! -f "$VCF_FILE" ]]; then
    echo "ERROR: VCF file not found: $VCF_FILE"
    exit 1
fi

# Backup old legend if it exists
if [[ -f "${LEGEND_OUT}.gz" ]]; then
    echo "Backing up old legend file..."
    mv "${LEGEND_OUT}.gz" "${LEGEND_OUT}.gz.backup"
    echo "✓ Backup created: ${LEGEND_OUT}.gz.backup"
    echo ""
fi

echo "Creating new legend file with proper variant IDs..."
echo "Input VCF: $VCF_FILE"
echo "Output: $LEGEND_OUT"
echo ""

# Create legend header
echo "id position a0 a1" > "$LEGEND_OUT"

# Extract variants and construct IDs in format: chr:pos:ref:alt
# Using colon separator to match chr1-22 format
zcat "$VCF_FILE" | grep -v "^#" | awk -v chr="$CHR" '{
    # Format: chr:position:ref:alt (using colons like chr1-22)
    id = chr ":" $2 ":" $4 ":" $5
    print id "\t" $2 "\t" $4 "\t" $5
}' >> "$LEGEND_OUT"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to create legend file"
    exit 1
fi

# Compress the legend file
echo "Compressing legend file..."
gzip "$LEGEND_OUT"

echo ""
echo "=========================================="
echo "✓ Legend file fixed successfully!"
echo "=========================================="
echo ""

# Show statistics
LEGEND_GZ="${LEGEND_OUT}.gz"
echo "Output file: $LEGEND_GZ"
echo "File size: $(du -h "$LEGEND_GZ" | cut -f1)"
echo "Variant count: $(zcat "$LEGEND_GZ" | wc -l | awk '{print $1-1}')"
echo ""

# Show first few lines
echo "First 10 variants:"
zcat "$LEGEND_GZ" | head -11
echo ""

echo "Next steps:"
echo "1. Regenerate hapbin files: sbatch scripts/ThisWorks/convert_fixed_chrX_to_hapbin.sh"
echo "2. Delete old iHS output: rm results/ihs/ALL.chrX.ihs.tsv"
echo "3. Rerun iHS: sbatch --array=22 scripts/ThisWorks/ihs_all.slurm"
