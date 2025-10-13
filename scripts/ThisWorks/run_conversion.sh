#!/usr/bin/env bash
# run_conversion.sh - Wrapper to convert IMPUTE to hapbin format

DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="$DATA_ROOT/hapbin"

mkdir -p "$OUT_DIR"

for CHR in "$@"; do
    echo "=========================================="
    echo "Processing chromosome $CHR..."
    echo "=========================================="
    
    # Find HAP file in multiple possible locations
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
        echo ""
        continue
    fi
    
    ANC="$DATA_ROOT/raw/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${CHR}.fa"
    MAP="$DATA_ROOT/maps/chr${CHR}.decode.sexavg.cm.tsv"
    
    # Check if ancestral FASTA exists
    if [[ ! -f "$ANC" ]]; then
        echo "  ✗ Ancestral FASTA not found: $ANC"
        echo ""
        continue
    fi
    
    # Auto-detect the exact chromosome name from the FASTA file
    # Extract chromosome name from the FASTA index or header
    if [[ -f "$ANC.fai" ]]; then
        FASTA_CHR=$(head -1 "$ANC.fai" | awk '{print $1}')
    else
        # Fallback: extract from FASTA header line
        FASTA_CHR=$(grep "^>" "$ANC" | head -1 | sed 's/^>//' | awk '{print $1}')
    fi
    
    echo "Using FASTA chromosome name: $FASTA_CHR"
    echo ""
    
    # Run the conversion
    python3 "$SCRIPT_DIR/convert_to_hapbin.py" \
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
        # Show output file sizes
        if [[ -f "$OUT_DIR/chr${CHR}.hapbin.hap" ]]; then
            HAP_SIZE=$(du -h "$OUT_DIR/chr${CHR}.hapbin.hap" | cut -f1)
            HAP_LINES=$(wc -l < "$OUT_DIR/chr${CHR}.hapbin.hap")
            echo "  Output: $HAP_LINES SNPs ($HAP_SIZE)"
        fi
    else
        echo ""
        echo "  ✗ Failed to convert chr$CHR"
    fi
    echo ""
done

echo "=========================================="
echo "CONVERSION COMPLETE"
echo "=========================================="