#!/usr/bin/env bash
#SBATCH --job-name=chrX_orig_to_hapbin
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/chrX_orig_to_hapbin_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/chrX_orig_to_hapbin_%j.err
#SBATCH --time=6:00:00
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
echo "Converting ORIGINAL chrX HAP to hapbin"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================="
echo ""
echo "The original chrX.impute.hap.gz ALREADY has:"
echo "  - 6404 fields on every line (diploid for all)"
echo "  - This is what hapbin needs!"
echo ""
echo "We just need to convert it to hapbin binary format."
echo ""

# Use the ORIGINAL HAP file (not the .fixed one!)
HAP_IN="$WORK_DIR/chr${CHR}.impute.hap.gz"
LEG_IN="$WORK_DIR/chr${CHR}.impute.legend.gz"

# Check inputs
echo "Checking input files..."
for file in "$HAP_IN" "$LEG_IN"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Missing input file: $file"
        exit 1
    fi
    echo "  ✓ $(basename $file)"
done
echo ""

# Verify the HAP file has consistent field counts (just check first line)
echo "Verifying HAP file format..."
# Temporarily disable pipefail for this command to avoid SIGPIPE issues
set +o pipefail
FIELDS_LINE1=$(zcat "$HAP_IN" | head -1 | awk '{print NF}')
set -o pipefail

echo "  Line 1: $FIELDS_LINE1 fields"

if [[ "$FIELDS_LINE1" != "6404" ]]; then
    echo ""
    echo "ERROR: Input HAP file does not have 6404 fields on first line!"
    echo "Expected 6404, got $FIELDS_LINE1"
    exit 1
fi

echo "  ✓ First line has 6404 fields (assuming consistent throughout)"
echo ""

# Prepare for conversion
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
    BACKUP_SUFFIX="backup_$(date +%Y%m%d_%H%M%S)"
    echo "Backing up old hapbin files..."
    mv "$OUT_DIR/chr${CHR}.hapbin.hap" "$OUT_DIR/chr${CHR}.hapbin.hap.$BACKUP_SUFFIX"
    mv "$OUT_DIR/chr${CHR}.hapbin.map" "$OUT_DIR/chr${CHR}.hapbin.map.$BACKUP_SUFFIX"
    echo "  Old files backed up with suffix: $BACKUP_SUFFIX"
    echo ""
fi

# Convert ORIGINAL HAP to hapbin format
echo "=========================================="
echo "Converting to hapbin format..."
echo "=========================================="
echo ""

# Run conversion using pixi from the project directory
PROJECT_DIR="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"

(cd "$PROJECT_DIR" && pixi run python "$SCRIPT_DIR/convert_to_hapbin.py" \
    --hap "$HAP_IN" \
    --legend "$LEG_IN" \
    --anc-fasta "$ANC" \
    --chr "$FASTA_CHR" \
    --recomb-tsv "$MAP" \
    --no-header \
    --pos-col 1 \
    --cm-col 2 \
    --sep ' ' \
    --out-prefix "$OUT_DIR/chr${CHR}.hapbin" \
    --verbose)

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to convert to hapbin"
    exit 1
fi

echo ""
echo "=========================================="
echo "✓ chrX conversion complete!"
echo "=========================================="
echo ""
echo "Output files:"
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
    echo "Checking output format:"
    head -1 "$OUT_DIR/chr${CHR}.hapbin.hap" | awk '{print "  Line 1:", NF, "fields"}'
    echo ""
fi

echo "Next step: rerun iHS calculation"
echo ""
echo "Command:"
echo "  pixi run ihsbin --hap $OUT_DIR/chr${CHR}.hapbin.hap \\"
echo "                  --map $OUT_DIR/chr${CHR}.hapbin.map \\"
echo "                  --out results/ihs/ALL.chr${CHR}.ihs.tsv"
echo ""
echo "Expected: iHS output should include PAR1 variants!"
