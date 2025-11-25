#!/usr/bin/env bash
#SBATCH --job-name=fix_chrX_hapbin_v2
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/fix_chrX_hapbin_v2_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/fix_chrX_hapbin_v2_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal

set -euo pipefail

CHR="X"
DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"
SCRIPT_DIR="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks"
WORK_DIR="$DATA_ROOT/work/chr${CHR}"
OUT_DIR="$DATA_ROOT/hapbin_with_dummy"

echo "=========================================="
echo "Fixing chromosome X for hapbin compatibility"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================="
echo ""
echo "This version creates HAP files with:"
echo "  - Constant field count (6404) on all lines"
echo "  - Males in non-PAR: real haplotype + '-' (dummy)"
echo "  - Females: haplotype1 + haplotype2"
echo "  - PAR regions: everyone gets 2 real haplotypes"
echo ""

# Input files
VCF="$DATA_ROOT/raw/chr${CHR}/chr${CHR}.vcf.gz"
HAP_IN="$WORK_DIR/chr${CHR}.impute.hap.gz"
LEG_IN="$WORK_DIR/chr${CHR}.impute.legend.gz"

# Fixed output
HAP_FIXED_HAPBIN="$WORK_DIR/chr${CHR}.impute.hap.hapbin_format.gz"

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

# Run the NEW fix script (for hapbin compatibility)
echo "Running chrX haplotype fix for hapbin..."
echo "This will:"
echo "  1. Infer sex from VCF genotypes"
echo "  2. Create HAP file with constant field count (6404)"
echo "  3. Use dummy '-' for male haplotypes in non-PAR"
echo ""

python3 "$SCRIPT_DIR/fix_chrX_haplotypes_for_hapbin.py" \
    --vcf "$VCF" \
    --hap-in "$HAP_IN" \
    --legend-in "$LEG_IN" \
    --hap-out "${HAP_FIXED_HAPBIN%.gz}" \
    --verbose

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to fix haplotypes"
    exit 1
fi

# Compress output
echo ""
echo "Compressing fixed HAP file..."
gzip -f "${HAP_FIXED_HAPBIN%.gz}"

# Verify output
echo ""
echo "Verifying output format..."
echo "Checking first line (PAR1):"
zcat "$HAP_FIXED_HAPBIN" | head -1 | awk '{print "  Fields:", NF}'

echo "Checking line 100000 (non-PAR):"
zcat "$HAP_FIXED_HAPBIN" | sed -n '100000p' | awk '{print "  Fields:", NF}'

echo ""
echo "Counting '-' characters in non-PAR line:"
zcat "$HAP_FIXED_HAPBIN" | sed -n '100000p' | tr ' ' '\n' | grep -c '^-$' || echo "  (count command failed, but may be OK)"

# Show statistics
echo ""
echo "Output file statistics:"
echo "  Size: $(du -h "$HAP_FIXED_HAPBIN" | cut -f1)"
echo "  Lines: $(zcat "$HAP_FIXED_HAPBIN" | wc -l)"
echo ""

# Now convert to hapbin format
echo "=========================================="
echo "Converting to hapbin format"
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

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"
echo "Output directory: $OUT_DIR"
echo ""

# Run conversion with FIXED haplotypes (hapbin-compatible format)
python3 "$SCRIPT_DIR/convert_to_hapbin.py" \
    --hap "$HAP_FIXED_HAPBIN" \
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
echo "  Fixed HAP (hapbin format): $HAP_FIXED_HAPBIN"
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
    echo "Verifying consistent field count:"
    head -1 "$OUT_DIR/chr${CHR}.hapbin.hap" | awk '{print "  Line 1 (PAR1):", NF, "fields"}'
    sed -n '100000p' "$OUT_DIR/chr${CHR}.hapbin.hap" | awk '{print "  Line 100000 (non-PAR):", NF, "fields"}'
    echo ""
    echo "Field counts should be identical (6404) for all lines"
fi

echo ""
echo "Next step: rerun iHS calculation with MODIFIED ihsbin that handles dummy '-' characters"
echo "  The iHS output should now include PAR1 regions!"
echo ""
echo "Command:"
echo "  /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/hapbin_source/build_custom/ihsbin \\"
echo "    --hap $OUT_DIR/chr${CHR}.hapbin.hap \\"
echo "    --map $OUT_DIR/chr${CHR}.hapbin.map \\"
echo "    --out results/ihs_with_dummy/ALL.chr${CHR}.ihs.tsv"
