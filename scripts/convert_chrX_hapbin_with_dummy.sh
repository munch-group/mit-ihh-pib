#!/usr/bin/env bash
#SBATCH --job-name=convert_chrX_hapbin
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/convert_chrX_hapbin_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/convert_chrX_hapbin_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal

set -euo pipefail

cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

CHR="X"
DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"
SCRIPT_DIR="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/ThisWorks"
WORK_DIR="$DATA_ROOT/work/chr${CHR}"
OUT_DIR="$DATA_ROOT/hapbin_with_dummy"

HAP_FIXED_HAPBIN="$WORK_DIR/chr${CHR}.impute.hap.hapbin_format.gz"
LEG_IN="$WORK_DIR/chr${CHR}.impute.legend.gz"
ANC="$DATA_ROOT/raw/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${CHR}.fa"
MAP="$DATA_ROOT/maps/chr${CHR}.decode.sexavg.cm.tsv"

echo "Converting chrX HAP with dummy variables to hapbin format..."
echo "Input: $HAP_FIXED_HAPBIN"
echo "Output: $OUT_DIR/"
echo ""

# Auto-detect chromosome name in FASTA
if [[ -f "$ANC.fai" ]]; then
    FASTA_CHR=$(head -1 "$ANC.fai" | awk '{print $1}')
else
    FASTA_CHR=$(grep "^>" "$ANC" | head -1 | sed 's/^>//' | awk '{print $1}')
fi

echo "Using FASTA chromosome name: $FASTA_CHR"
echo ""

# Create output directory
mkdir -p "$OUT_DIR"

# Run conversion with FIXED haplotypes (hapbin-compatible format with dummy '-')
pixi run python "$SCRIPT_DIR/convert_to_hapbin.py" \
    --hap "$HAP_FIXED_HAPBIN" \
    --legend "$LEG_IN" \
    --anc-fasta "$ANC" \
    --chr "$FASTA_CHR" \
    --recomb-tsv "$MAP" \
    --no-header \
    --pos-col 1 \
    --cm-col 2 \
    --sep ' ' \
    --out-prefix "$OUT_DIR/chr${CHR}.hapbin" \
    --verbose

echo ""
echo "✓ Conversion complete!"
echo "Output files:"
echo "  HAP: $OUT_DIR/chr${CHR}.hapbin.hap"
echo "  MAP: $OUT_DIR/chr${CHR}.hapbin.map"
