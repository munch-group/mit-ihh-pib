#!/usr/bin/env bash

# This script performs both the cumulative cM calculation (with 2/3 adjustment)
# AND the final reformatting step required by Relate's --map flag.

# --- Define Paths ---
DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"
MAPS_DIR="$DATA_ROOT/maps"

BIGWIG_FILE="$MAPS_DIR/recombAvg.hg38.bw"
BIGWIG_TO_BEDGRAPH="$MAPS_DIR/bigWigToBedGraph"

BEDGRAPH_FILE="$MAPS_DIR/chrX.decode.sexavg.rate.bedGraph"
CUMULATIVE_CM_FILE="$MAPS_DIR/chrX.decode.sexavg.cm.tsv"
OUTPUT_RELATE_MAP="$MAPS_DIR/chrX.relate.map" # The final file for Relate

echo "=========================================="
echo "Creating Relate-Compatible Map for chrX"
echo "=========================================="

# Check if necessary tools/files exist
if [[ ! -f "$BIGWIG_FILE" ]]; then
    echo "ERROR: BigWig file not found: $BIGWIG_FILE"
    exit 1
fi

if [[ ! -f "$BIGWIG_TO_BEDGRAPH" ]]; then
    echo "ERROR: bigWigToBedGraph tool not found: $BIGWIG_TO_BEDGRAPH"
    exit 1
fi

# ----------------------------------------------------
# STEP 1: Extract BedGraph and Calculate Cumulative cM
# ----------------------------------------------------

if [[ ! -f "$CUMULATIVE_CM_FILE" ]]; then
    echo "Extracting chrX from BigWig to BedGraph and calculating cumulative cM..."
    
    # Extract BedGraph, or use existing one
    if [[ ! -f "$BEDGRAPH_FILE" ]]; then
        echo "   -> Running bigWigToBedGraph..."
        "$BIGWIG_TO_BEDGRAPH" "$BIGWIG_FILE" "$BEDGRAPH_FILE" -chrom=chrX
        if [[ $? -ne 0 ]]; then
            echo "ERROR: Failed to extract chrX to BedGraph"
            exit 1
        fi
    else
        echo "   -> BedGraph file already exists: $BEDGRAPH_FILE (Skipping extraction)"
    fi
    
    echo "   -> Converting BedGraph to cumulative cM (2/3 adjusted)..."
    
    # AWK block for cumulative cM calculation
    awk 'BEGIN { cumulative_cm = 0 }
    {
        start = $2
        end = $3
        rate_per_mb = $4  # cM per Mb
        
        distance_mb = (end - start) / 1000000.0
        cm_interval = rate_per_mb * distance_mb
        cumulative_cm += cm_interval
        
        # Apply 2/3 adjustment for X chromosome
        adjusted_cm = cumulative_cm * 2.0/3.0
        
        # Output: chr, position (end of bin), cM
        print $1, end, adjusted_cm
    }' "$BEDGRAPH_FILE" > "$CUMULATIVE_CM_FILE"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to create cumulative cM file"
        exit 1
    fi
    echo "✓ Created cumulative cM file: $CUMULATIVE_CM_FILE"
else
    echo "Cumulative cM file already exists: $CUMULATIVE_CM_FILE (Skipping Step 1)"
fi


# ----------------------------------------------------
# STEP 2: Reformat to Relate's 3-column Standard
# ----------------------------------------------------

echo ""
echo "Reformatting to Relate's 3-column, space-separated standard..."

# The most robust single-line AWK command to avoid all parsing errors.
# Input is CUMULATIVE_CM_FILE (Col 1: Chr, Col 2: Pos, Col 3: cM)
# Output is RELATE_MAP_FILE (Col 1: ID, Col 2: Pos, Col 3: cM)
awk '{print $1 ":" $2, $2, $3}' "$CUMULATIVE_CM_FILE" > "$OUTPUT_RELATE_MAP"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to create Relate map file: $OUTPUT_RELATE_MAP"
    exit 1
fi

echo "✓ Created Relate map file: $OUTPUT_RELATE_MAP"

# ----------------------------------------------------
# Final Statistics
# ----------------------------------------------------
echo ""
echo "=========================================="
echo "Final Statistics"
echo "=========================================="
LINE_COUNT=$(wc -l < "$OUTPUT_RELATE_MAP")
FINAL_CM=$(tail -1 "$OUTPUT_RELATE_MAP" | awk '{print $3}')

echo "  Total positions: $LINE_COUNT"
echo "  Final genetic length: $FINAL_CM cM"
echo ""
echo "First 5 lines of the final Relate Map (chrX.relate.map):"
head -5 "$OUTPUT_RELATE_MAP"
echo ""
echo "You can now use this file with the Relate command:"
echo "  ... --map $OUTPUT_RELATE_MAP ..."