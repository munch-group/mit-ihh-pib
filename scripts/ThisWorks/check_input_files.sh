#!/usr/bin/env bash
# Check which chromosomes have all necessary input files

DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"

echo "Checking available files for chromosomes 1-22..."
echo ""
echo "Chr | HAP/LEGEND | Ancestral FASTA | Recomb Map | Status"
echo "----+------------+-----------------+------------+--------"

for CHR in {1..22}; do
    HAS_HAP=0
    HAS_ANC=0
    HAS_MAP=0
    
    # Check for HAP/LEGEND in multiple locations
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
    
    # Determine status
    if [[ $HAS_HAP -eq 1 && $HAS_ANC -eq 1 && $HAS_MAP -eq 1 ]]; then
        STATUS="✓ Ready"
    else
        STATUS="✗ Missing:"
        [[ $HAS_HAP -eq 0 ]] && STATUS="$STATUS HAP"
        [[ $HAS_ANC -eq 0 ]] && STATUS="$STATUS ANC"
        [[ $HAS_MAP -eq 0 ]] && STATUS="$STATUS MAP"
    fi
    
    printf "%3d | %-10s | %-15s | %-10s | %s\n" \
        $CHR \
        $([[ $HAS_HAP -eq 1 ]] && echo "✓" || echo "✗") \
        $([[ $HAS_ANC -eq 1 ]] && echo "✓" || echo "✗") \
        $([[ $HAS_MAP -eq 1 ]] && echo "✓" || echo "✗") \
        "$STATUS"
done

echo ""
echo "Already converted:"
if [[ -d "$DATA_ROOT/hapbin" ]]; then
    ls "$DATA_ROOT/hapbin/"chr*.hapbin.hap 2>/dev/null | while read f; do
        CHR=$(basename "$f" .hapbin.hap | sed 's/chr//')
        LINES=$(wc -l < "$f")
        SIZE=$(du -h "$f" | cut -f1)
        echo "  Chr$CHR: $LINES SNPs ($SIZE)"
    done
else
    echo "  None yet"
fi