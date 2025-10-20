#!/usr/bin/env bash
#SBATCH --job-name=fix_chrX_hapbin
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/fix_chrX_hapbin_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/fix_chrX_hapbin_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal

set -euo pipefail

CHR="X"
DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"
SCRIPT_DIR="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks"
WORK_DIR="$DATA_ROOT/work/chr${CHR}"
OUT_DIR="$DATA_ROOT/hapbin"

echo "=========================================="
echo "Fixing chromosome X haplotype data"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================="
echo ""

# Input files
VCF="$DATA_ROOT/raw/chr${CHR}/chr${CHR}.vcf.gz"
HAP_IN="$WORK_DIR/chr${CHR}.impute.hap.gz"
LEG_IN="$WORK_DIR/chr${CHR}.impute.legend.gz"

# Fixed output
HAP_FIXED="$WORK_DIR/chr${CHR}.impute.hap.fixed.gz"

# Check inputs
echo "Checking input files..."
for file in "$VCF" "$HAP_IN" "$LEG_IN"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Missing input file: $file"
        exit 1
    fi
    echo "  ✓ $(basename $file)"
done
echo ""

# Run the fix script
echo "Running chrX haplotype fix..."
echo "This will:"
echo "  1. Infer sex from VCF genotypes"
echo "  2. Remove duplicate male haplotypes in non-PAR regions"
echo "  3. Keep both haplotypes for females and PAR regions"
echo ""

python3 "$SCRIPT_DIR/fix_chrX_haplotypes.py" \
    --vcf "$VCF" \
    --hap-in "$HAP_IN" \
    --legend-in "$LEG_IN" \
    --hap-out "${HAP_FIXED%.gz}" \
    --verbose

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to fix haplotypes"
    exit 1
fi

# Compress output
echo ""
echo "Compressing fixed HAP file..."
gzip -f "${HAP_FIXED%.gz}"

# Show statistics
echo ""
echo "Checking output..."
echo "  Fixed HAP file: $(du -h "$HAP_FIXED" | cut -f1)"
echo "  Note: Haplotype count varies by region (PAR vs non-PAR)"
echo "  PAR regions: ~6404 haplotypes (diploid for all)"
echo "  Non-PAR regions: ~5968 haplotypes (haploid for males)"
echo ""

# Now convert to hapbin format
echo "=========================================="
echo "Converting fixed HAP to hapbin format"
echo "=========================================="
echo ""

ANC="$DATA_ROOT/raw/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${CHR}.fa"
MAP="$DATA_ROOT/maps/chr${CHR}.decode.sexavg.cm.tsv"

if [[ ! -f "$ANC" ]]; then
    echo "ERROR: Ancestral FASTA not found: $ANC"
    exit 1
fi

if [[ ! -f "$MAP" ]]; then
    echo "ERROR: Recombination map not found: $MAP"
    exit 1
fi

# Auto-detect chromosome name in FASTA
if [[ -f "$ANC.fai" ]]; then
    FASTA_CHR=$(head -1 "$ANC.fai" | awk '{print $1}')
else
    FASTA_CHR=$(grep "^>" "$ANC" | head -1 | sed 's/^>//' | awk '{print $1}')
fi

echo "Using FASTA chromosome name: $FASTA_CHR"
echo ""

# Backup old hapbin files if they exist
if [[ -f "$OUT_DIR/chr${CHR}.hapbin.hap" ]]; then
    echo "Backing up old hapbin files..."
    mv "$OUT_DIR/chr${CHR}.hapbin.hap" "$OUT_DIR/chr${CHR}.hapbin.hap.backup"
    mv "$OUT_DIR/chr${CHR}.hapbin.map" "$OUT_DIR/chr${CHR}.hapbin.map.backup"
    echo "  Old files backed up with .backup extension"
    echo ""
fi

# Run conversion with FIXED haplotypes
python3 "$SCRIPT_DIR/convert_to_hapbin.py" \
    --hap "$HAP_FIXED" \
    --legend "$LEG_IN" \
    --anc-fasta "$ANC" \
    --chr "$FASTA_CHR" \
    --recomb-tsv "$MAP" \
    --no-header \
    --pos-col 1 \
    --cm-col 2 \
    --out-prefix "$OUT_DIR/chr${CHR}.hapbin" \
    --verbose

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to convert to hapbin"
    exit 1
fi

echo ""
echo "=========================================="
echo "✓ chrX processing complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Fixed HAP: $HAP_FIXED"
echo "  hapbin HAP: $OUT_DIR/chr${CHR}.hapbin.hap"
echo "  hapbin MAP: $OUT_DIR/chr${CHR}.hapbin.map"
echo ""

if [[ -f "$OUT_DIR/chr${CHR}.hapbin.hap" ]]; then
    HAP_SIZE=$(du -h "$OUT_DIR/chr${CHR}.hapbin.hap" | cut -f1)
    HAP_LINES=$(wc -l < "$OUT_DIR/chr${CHR}.hapbin.hap")
    echo "Final hapbin file:"
    echo "  Size: $HAP_SIZE"
    echo "  SNPs: $HAP_LINES"
    echo ""
    echo "Note: Chromosome count varies by region (PAR vs non-PAR)"
    echo "  Was 6404 for all variants before fix"
    echo "  Now ~6404 in PAR, ~5968 in non-PAR"
fi

echo "Next step: rerun iHS calculation"
echo "  pixi run ihsbin --hap $OUT_DIR/chr${CHR}.hapbin.hap --map $OUT_DIR/chr${CHR}.hapbin.map --out results/ihs/ALL.chr${CHR}.ihs.tsv"
