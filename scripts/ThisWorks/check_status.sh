#!/usr/bin/env bash
# check_status.sh - Check what we have and what's missing

DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"
OUT_DIR="$DATA_ROOT/hapbin"

echo "=========================================="
echo "CHROMOSOME STATUS CHECK"
echo "=========================================="
echo ""
echo "Checking chromosomes 1-22..."
echo ""

# Arrays to track status
declare -a READY=()
declare -a MISSING_HAP=()
declare -a MISSING_ANC=()
declare -a MISSING_MAP=()
declare -a CONVERTED=()
declare -a FAILED=()

for CHR in {1..22}; do
    HAS_HAP=0
    HAS_ANC=0
    HAS_MAP=0
    HAS_OUTPUT=0
    
    # Check for HAP/LEGEND
    if [[ -f "$DATA_ROOT/work/chr${CHR}/chr${CHR}.impute.hap.gz" ]] || \
       [[ -f "$DATA_ROOT/raw/chr${CHR}/chr${CHR}.impute.hap.gz" ]] || \
       [[ -f "$DATA_ROOT/impute/chr${CHR}.hap.gz" ]]; then
        HAS_HAP=1
    fi
    
    # Check for ancestral FASTA
    if [[ -f "$DATA_ROOT/raw/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${CHR}.fa" ]]; then
        HAS_ANC=1
    fi
    
    # Check for recombination map
    if [[ -f "$DATA_ROOT/maps/chr${CHR}.decode.sexavg.cm.tsv" ]]; then
        HAS_MAP=1
    fi
    
    # Check for output files
    if [[ -f "$OUT_DIR/chr${CHR}.hapbin.hap" ]] && [[ -f "$OUT_DIR/chr${CHR}.hapbin.map" ]]; then
        HAS_OUTPUT=1
    fi
    
    # Categorize
    if [[ $HAS_OUTPUT -eq 1 ]]; then
        CONVERTED+=($CHR)
    elif [[ $HAS_HAP -eq 1 && $HAS_ANC -eq 1 && $HAS_MAP -eq 1 ]]; then
        READY+=($CHR)
    else
        if [[ $HAS_HAP -eq 0 ]]; then
            MISSING_HAP+=($CHR)
        fi
        if [[ $HAS_ANC -eq 0 ]]; then
            MISSING_ANC+=($CHR)
        fi
        if [[ $HAS_MAP -eq 0 ]]; then
            MISSING_MAP+=($CHR)
        fi
    fi
done

echo "✓ ALREADY CONVERTED (${#CONVERTED[@]}):"
if [[ ${#CONVERTED[@]} -gt 0 ]]; then
    for CHR in "${CONVERTED[@]}"; do
        HAP_SIZE=$(du -h "$OUT_DIR/chr${CHR}.hapbin.hap" 2>/dev/null | cut -f1)
        HAP_LINES=$(wc -l < "$OUT_DIR/chr${CHR}.hapbin.hap" 2>/dev/null)
        echo "  Chr${CHR}: $HAP_LINES SNPs ($HAP_SIZE)"
    done
else
    echo "  None"
fi
echo ""

echo "🟢 READY TO CONVERT (${#READY[@]}):"
if [[ ${#READY[@]} -gt 0 ]]; then
    echo "  Chromosomes: ${READY[*]}"
else
    echo "  None"
fi
echo ""

echo "🔴 MISSING INPUT FILES:"
echo ""

if [[ ${#MISSING_HAP[@]} -gt 0 ]]; then
    echo "  Missing HAP/LEGEND (${#MISSING_HAP[@]}):"
    echo "    Chromosomes: ${MISSING_HAP[*]}"
    echo ""
fi

if [[ ${#MISSING_ANC[@]} -gt 0 ]]; then
    echo "  Missing Ancestral FASTA (${#MISSING_ANC[@]}):"
    echo "    Chromosomes: ${MISSING_ANC[*]}"
    echo ""
fi

if [[ ${#MISSING_MAP[@]} -gt 0 ]]; then
    echo "  Missing Recombination Map (${#MISSING_MAP[@]}):"
    echo "    Chromosomes: ${MISSING_MAP[*]}"
    echo ""
fi

echo "=========================================="
echo "SLURM JOB LOG CHECK"
echo "=========================================="
echo ""

# Check for SLURM logs
LOGS_DIR="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs"

if [[ -d "$LOGS_DIR" ]]; then
    echo "Recent SLURM jobs:"
    ls -lht "$LOGS_DIR"/hapbin_chr* 2>/dev/null | head -20
    echo ""
    
    echo "Checking for errors in logs..."
    for CHR in {1..22}; do
        ERR_FILE=$(ls -t "$LOGS_DIR"/hapbin_chr${CHR}_*.err 2>/dev/null | head -1)
        if [[ -f "$ERR_FILE" ]]; then
            if [[ -s "$ERR_FILE" ]]; then
                echo "  Chr${CHR}: Errors found in $ERR_FILE"
                echo "    $(tail -3 "$ERR_FILE" | sed 's/^/      /')"
            fi
        fi
    done
else
    echo "No logs directory found at: $LOGS_DIR"
fi

echo ""
echo "=========================================="
echo "SUMMARY"
echo "=========================================="
echo "Converted:    ${#CONVERTED[@]}/22"
echo "Ready:        ${#READY[@]}/22"
echo "Missing data: $((${#MISSING_HAP[@]} + ${#MISSING_ANC[@]} + ${#MISSING_MAP[@]}))"
echo ""

if [[ ${#READY[@]} -gt 0 ]]; then
    echo "💡 To convert ready chromosomes, run:"
    echo "   sbatch --array=$(IFS=,; echo "${READY[*]}") scripts/ThisWorks/run_conversion_array.slurm"
    echo ""
fi

if [[ ${#MISSING_HAP[@]} -gt 0 ]] || [[ ${#MISSING_ANC[@]} -gt 0 ]] || [[ ${#MISSING_MAP[@]} -gt 0 ]]; then
    echo "⚠️  Need to download missing files for: ${MISSING_HAP[*]} ${MISSING_ANC[*]} ${MISSING_MAP[*]}"
fi