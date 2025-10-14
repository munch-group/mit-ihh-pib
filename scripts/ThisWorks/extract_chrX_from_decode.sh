#!/usr/bin/env bash

# Extract chrX from BigWig and create .cm.tsv with 2/3 adjustment
# Following the same process used for chr1-22

DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"
MAPS_DIR="$DATA_ROOT/maps"

BIGWIG_FILE="$MAPS_DIR/recombAvg.hg38.bw"
BIGWIG_TO_BEDGRAPH="$MAPS_DIR/bigWigToBedGraph"

BEDGRAPH_FILE="$MAPS_DIR/chrX.decode.sexavg.rate.bedGraph"
OUTPUT_MAP="$MAPS_DIR/chrX.decode.sexavg.cm.tsv"

echo "=========================================="
echo "Creating chrX recombination map"
echo "=========================================="

# Check if files exist
if [[ ! -f "$BIGWIG_FILE" ]]; then
    echo "ERROR: BigWig file not found: $BIGWIG_FILE"
    exit 1
fi

if [[ ! -f "$BIGWIG_TO_BEDGRAPH" ]]; then
    echo "ERROR: bigWigToBedGraph tool not found: $BIGWIG_TO_BEDGRAPH"
    exit 1
fi

# Step 1: Extract chrX to BedGraph (same as for other chromosomes)
if [[ ! -f "$BEDGRAPH_FILE" ]]; then
    echo "Extracting chrX from BigWig to BedGraph..."
    "$BIGWIG_TO_BEDGRAPH" "$BIGWIG_FILE" "$BEDGRAPH_FILE" -chrom=chrX
    
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to extract chrX"
        exit 1
    fi
    echo "✓ Created: $BEDGRAPH_FILE"
else
    echo "BedGraph file already exists: $BEDGRAPH_FILE"
fi

# Step 2: Convert BedGraph to cumulative cM format with 2/3 adjustment
echo ""
echo "Converting BedGraph to cumulative cM format..."
echo "Applying 2/3 adjustment for X chromosome..."

awk 'BEGIN {
    cumulative_cm = 0
}
{
    chr = $1
    start = $2
    end = $3
    rate_per_mb = $4  # cM per Mb
    
    # Calculate distance in Mb
    distance_mb = (end - start) / 1000000.0
    
    # Calculate cM for this interval
    cm_interval = rate_per_mb * distance_mb
    
    # Add to cumulative
    cumulative_cm += cm_interval
    
    # Apply 2/3 adjustment for X chromosome
    # (because only females recombine on X, not males)
    adjusted_cm = cumulative_cm * 2.0/3.0
    
    # Output in same format as other chromosomes: chr position cM
    print chr, end, adjusted_cm
}' "$BEDGRAPH_FILE" > "$OUTPUT_MAP"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to create .cm.tsv file"
    exit 1
fi

echo "✓ Created: $OUTPUT_MAP"
echo ""

# Show statistics
LINE_COUNT=$(wc -l < "$OUTPUT_MAP")
FINAL_CM=$(tail -1 "$OUTPUT_MAP" | awk '{print $3}')

echo "Statistics:"
echo "  Total positions: $LINE_COUNT"
echo "  Total genetic length (2/3 adjusted): $FINAL_CM cM"
echo ""
echo "First 10 lines:"
head -10 "$OUTPUT_MAP"
echo ""
echo "Last 5 lines:"
tail -5 "$OUTPUT_MAP"

echo ""
echo "=========================================="
echo "✓ Complete!"
echo "=========================================="
echo ""
echo "The chrX map is now ready in the same format as chr1-22:"
echo "  $OUTPUT_MAP"
echo ""
echo "Next step: Run vcf2hap_chrX.sh to process the VCF data"