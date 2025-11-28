#!/usr/bin/env bash
# Reformat iHS output to match the older format
# Removes Index and Freq columns, renames ID to Location

set -euo pipefail

INPUT="$1"
OUTPUT="$2"

echo "Reformatting iHS output..."
echo "  Input:  $INPUT"
echo "  Output: $OUTPUT"
echo ""

# Extract columns: ID (col 2), iHH_0 (col 4), iHH_1 (col 5), iHS (col 6), Std iHS (col 7)
# Rename ID to Location in the header
awk 'BEGIN {FS="\t"; OFS="\t"}
     NR==1 {print "Location", $4, $5, $6, $7}
     NR>1 {print $2, $4, $5, $6, $7}' "$INPUT" > "$OUTPUT"

echo "✓ Reformatting complete"
echo ""
echo "Original columns: $(head -1 "$INPUT")"
echo "New columns:      $(head -1 "$OUTPUT")"
echo ""
echo "Line counts:"
echo "  Input:  $(wc -l < "$INPUT")"
echo "  Output: $(wc -l < "$OUTPUT")"
