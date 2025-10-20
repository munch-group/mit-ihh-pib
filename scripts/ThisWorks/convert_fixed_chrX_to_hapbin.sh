#!/usr/bin/env bash
#SBATCH --job-name=chrX_to_hapbin
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/chrX_to_hapbin_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/chrX_to_hapbin_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal

set -euo pipefail

CHR="X"
DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"
SCRIPT_DIR="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks"
WORK_DIR="$DATA_ROOT/work/chr${CHR}"
OUT_DIR="$DATA_ROOT/hapbin"

echo "=========================================="
echo "Converting fixed chrX HAP to hapbin format"
echo "=========================================="
echo ""

# Fixed HAP file (already exists from previous job)
HAP_FIXED="$WORK_DIR/chr${CHR}.impute.hap.fixed.gz"
LEG_IN="$WORK_DIR/chr${CHR}.impute.legend.gz"
ANC="$DATA_ROOT/raw/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${CHR}.fa"
MAP="$DATA_ROOT/maps/chr${CHR}.decode.sexavg.cm.tsv"

# Check inputs
if [[ ! -f "$HAP_FIXED" ]]; then
    echo "ERROR: Fixed HAP file not found: $HAP_FIXED"
    exit 1
fi

echo "Input files:"
echo "  Fixed HAP: $HAP_FIXED ($(du -h "$HAP_FIXED" | cut -f1))"
echo "  Legend: $LEG_IN"
echo "  Ancestral: $ANC"
echo "  Recomb map: $MAP"
echo ""

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
    echo "✓ Old files backed up"
    echo ""
fi

# Run conversion
echo "Running conversion to hapbin format..."
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
    echo "ERROR: Conversion failed"
    exit 1
fi

echo ""
echo "=========================================="
echo "✓ Conversion complete!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  $OUT_DIR/chr${CHR}.hapbin.hap"
echo "  $OUT_DIR/chr${CHR}.hapbin.map"
echo ""

if [[ -f "$OUT_DIR/chr${CHR}.hapbin.hap" ]]; then
    HAP_SIZE=$(du -h "$OUT_DIR/chr${CHR}.hapbin.hap" | cut -f1)
    echo "  Size: $HAP_SIZE"
    echo ""
fi

echo "Next: Delete old iHS results and rerun"
echo "  rm results/ihs/ALL.chrX.ihs.tsv"
echo "  sbatch --array=22 scripts/ThisWorks/ihs_all.slurm"
