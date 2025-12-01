#!/usr/bin/env bash
#SBATCH -J ihs-chrX-impute2
#SBATCH -A mit-ihh-pib
#SBATCH -p normal
#SBATCH -c 8
#SBATCH -t 36:00:00
#SBATCH -o /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_chrX_impute2_%j.out
#SBATCH -e /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/ihs_chrX_impute2_%j.err

set -euo pipefail

REPO_ROOT="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib"
DATA_DIR="/home/vanbruggenmit/mit-ihh-pib/data/grch38/work/chrX_impute2_conversion"

# Use modified ihsbin that supports '-' and '?' characters
IHSBIN="$REPO_ROOT/hapbin_source/build_custom/ihsbin"

cd "$REPO_ROOT"

CHR="X"

OUT_DIR="$REPO_ROOT/results/ihs_impute2"
mkdir -p "$OUT_DIR" "$REPO_ROOT/logs"

HAP="$DATA_DIR/chr${CHR}.impute.hap.gz"
MAP="$DATA_DIR/chr${CHR}.impute.map"
OUT="$OUT_DIR/ALL.chr${CHR}.ihs.tsv"

echo "=========================================="
echo "iHS calculation for chromosome ${CHR}"
echo "Using IMPUTE2 format files (gold standard)"
echo "Job ID: $SLURM_JOB_ID"
echo "=========================================="
echo "Inputs:"
echo "  HAP: $HAP"
echo "  MAP: $MAP"
echo "  IHSBIN: $IHSBIN"
echo "Output:"
echo "  OUT: $OUT"
echo ""

# Verify modified ihsbin exists
if [[ ! -x "$IHSBIN" ]]; then
    echo "ERROR: Modified ihsbin binary not found or not executable: $IHSBIN"
    echo "Please rebuild hapbin with: cd $REPO_ROOT/hapbin_source && ./build.sh"
    exit 1
fi
echo "✓ Modified ihsbin found and executable"
echo ""

# Check HAP file format
echo "Verifying HAP file format..."
HAP_LINES=$(zcat "$HAP" | wc -l)
HAP_SIZE=$(du -h "$HAP" | cut -f1)
FIELDS_LINE1=$(zcat "$HAP" | head -1 | awk '{print NF}')
echo "  HAP file: $HAP_SIZE, $HAP_LINES SNPs"
echo "  First line: $FIELDS_LINE1 fields (should be 6404)"
echo ""

# Check MAP file
echo "Verifying MAP file..."
MAP_LINES=$(wc -l < "$MAP")
MAP_FIRST_POS=$(head -1 "$MAP" | awk '{print $4}')
MAP_LAST_POS=$(tail -1 "$MAP" | awk '{print $4}')
echo "  MAP file: $MAP_LINES positions"
echo "  Position range: $MAP_FIRST_POS - $MAP_LAST_POS"
echo ""

# Validate inputs
if [[ ! -s "$HAP" ]]; then
    echo "ERROR: HAP file missing or empty: $HAP"
    exit 1
fi

if [[ ! -s "$MAP" ]]; then
    echo "ERROR: MAP file missing or empty: $MAP"
    exit 1
fi

if [[ "$HAP_LINES" != "$MAP_LINES" ]]; then
    echo "WARNING: HAP and MAP have different line counts"
    echo "  HAP: $HAP_LINES"
    echo "  MAP: $MAP_LINES"
fi

echo "Running iHS for chr${CHR}..."
echo "Parameters:"
echo "  --minmaf 0.05 (minimum minor allele frequency)"
echo "  --cutoff 0.05 (EHH decay cutoff)"
echo ""
echo "Note: Using modified ihsbin that handles '-' and '?' characters"
echo "      '-' = males' missing second haplotype in non-PAR region"
echo "      '?' = unknown/missing genotype data (from VCF missing genotypes)"
echo ""

# Run ihsbin with IMPUTE2 format files
# The modified ihsbin handles both '-' and '?' characters
"$IHSBIN" --hap "$HAP" --map "$MAP" --out "$OUT" --minmaf 0.05 --cutoff 0.05

echo ""
echo "=========================================="
if [[ -s "$OUT" ]]; then
    OUT_LINES=$(wc -l < "$OUT")
    OUT_SIZE=$(du -h "$OUT" | cut -f1)
    echo "✓ iHS completed for chr${CHR}"
    echo "=========================================="
    echo ""
    echo "Output file: $OUT"
    echo "  Size: $OUT_SIZE"
    echo "  Lines: $OUT_LINES"
    echo ""

    # Show first few results
    echo "First 5 iHS results:"
    head -6 "$OUT"
    echo ""

    # Check if PAR1 is included
    echo "Checking PAR1 coverage (expect variants starting ~95,187):"
    tail -n +2 "$OUT" | head -5 | awk '{print "  Position:", $2}'
    echo ""

    # Summary statistics
    echo "Summary statistics:"
    tail -n +2 "$OUT" | awk '
        BEGIN {count=0; sum_ihs=0; sum_abs_ihs=0}
        {
            count++
            sum_ihs += $5
            sum_abs_ihs += ($5 < 0 ? -$5 : $5)
        }
        END {
            print "  Valid loci:", count
            print "  Mean iHS:", sum_ihs/count
            print "  Mean |iHS|:", sum_abs_ihs/count
        }
    '
    echo ""
else
    echo "ERROR: Output file not created or is empty"
    echo "=========================================="
    exit 1
fi

echo "Expected: iHS output should include all PAR1, non-PAR, and PAR2 variants!"
echo "Total variants: 2,858,131"
echo ""