#!/usr/bin/env bash
# check_naming.sh - Diagnostic script to check file naming conventions

DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"

echo "=========================================="
echo "FILE NAMING DIAGNOSTIC"
echo "=========================================="
echo ""

for CHR in 21 22; do
    echo "==================== CHROMOSOME $CHR ===================="
    echo ""
    
    # 1. Check HAP/LEGEND files
    echo "1. HAP/LEGEND FILES:"
    if [[ -f "$DATA_ROOT/work/chr${CHR}/chr${CHR}.impute.hap.gz" ]]; then
        HAP="$DATA_ROOT/work/chr${CHR}/chr${CHR}.impute.hap.gz"
        LEG="$DATA_ROOT/work/chr${CHR}/chr${CHR}.impute.legend.gz"
        echo "   Found in: work/chr${CHR}/"
    elif [[ -f "$DATA_ROOT/raw/chr${CHR}/chr${CHR}.impute.hap.gz" ]]; then
        HAP="$DATA_ROOT/raw/chr${CHR}/chr${CHR}.impute.hap.gz"
        LEG="$DATA_ROOT/raw/chr${CHR}/chr${CHR}.impute.legend.gz"
        echo "   Found in: raw/chr${CHR}/"
    elif [[ -f "$DATA_ROOT/impute/chr${CHR}.hap.gz" ]]; then
        HAP="$DATA_ROOT/impute/chr${CHR}.hap.gz"
        LEG="$DATA_ROOT/impute/chr${CHR}.legend.gz"
        echo "   Found in: impute/"
    else
        echo "   ✗ NOT FOUND"
        continue
    fi
    echo "   HAP: $HAP"
    echo "   LEGEND: $LEG"
    echo ""
    
    # 2. Check LEGEND header
    echo "2. LEGEND HEADER:"
    if [[ -f "$LEG" ]]; then
        echo "   First line:"
        zcat "$LEG" 2>/dev/null | head -1
        echo ""
        echo "   First data line:"
        zcat "$LEG" 2>/dev/null | head -2 | tail -1
    fi
    echo ""
    
    # 3. Check ancestral FASTA
    echo "3. ANCESTRAL FASTA:"
    ANC="$DATA_ROOT/raw/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${CHR}.fa"
    if [[ -f "$ANC" ]]; then
        echo "   File: $ANC"
        echo "   Chromosome names in FASTA:"
        if command -v samtools &> /dev/null; then
            samtools faidx "$ANC" 2>/dev/null
            cat "$ANC.fai" 2>/dev/null | awk '{print "     - " $1 " (length: " $2 " bp)"}'
        else
            # Fallback: just grep for the header
            grep "^>" "$ANC" | head -5 | sed 's/^/     /'
        fi
    else
        echo "   ✗ NOT FOUND: $ANC"
    fi
    echo ""
    
    # 4. Check recombination map
    echo "4. RECOMBINATION MAP:"
    MAP="$DATA_ROOT/maps/chr${CHR}.decode.sexavg.cm.tsv"
    if [[ -f "$MAP" ]]; then
        echo "   File: $MAP"
        echo "   First 3 lines:"
        head -3 "$MAP" | nl -w2 -s': '
        echo ""
        echo "   Column structure:"
        head -1 "$MAP" | awk -F'\t' '{for(i=1;i<=NF;i++) print "     Col " i-1 ": " $i}'
    else
        echo "   ✗ NOT FOUND: $MAP"
    fi
    echo ""
    echo ""
done

echo "=========================================="
echo "SUMMARY & RECOMMENDATIONS"
echo "=========================================="
echo ""
echo "Based on the output above, we need to determine:"
echo "1. Which column in the recomb map is position (likely col 1)"
echo "2. Which column in the recomb map is cM (likely col 2)"
echo "3. Does the recomb map have a header row? (check if line 1 looks like data)"
echo "4. What is the exact chromosome name in the ancestral FASTA?"
echo ""