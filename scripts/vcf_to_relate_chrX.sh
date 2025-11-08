#!/usr/bin/env bash
#SBATCH --job-name=vcf_to_relate_chrX
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/vcf_to_relate_chrX_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/vcf_to_relate_chrX_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# Convert chromosome X VCF to Relate format
# Note: X chromosome requires special handling due to hemizygosity in males

# Set up Relate path
RELATE_BIN="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/software/relate/bin"
export PATH="$RELATE_BIN:$PATH"

# Check if Relate is available
if ! command -v RelateFileFormats &> /dev/null; then
    echo "ERROR: RelateFileFormats not found"
    echo "Please run: sbatch scripts/install_relate.sh"
    exit 1
fi

CHR="X"

echo "=========================================="
echo "Converting chromosome X to Relate format"
echo "=========================================="

# Paths
VCF_DIR="/faststorage/project/mit-ihh-pib/data/grch38/raw/chrX"
VCF_FILE="$VCF_DIR/chrX.vcf.gz"
OUTPUT_DIR="/faststorage/project/mit-ihh-pib/data/grch38/relate"
OUTPUT_PREFIX="$OUTPUT_DIR/chrX"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check if VCF exists
if [[ ! -f "$VCF_FILE" ]]; then
    echo "ERROR: VCF file not found: $VCF_FILE"
    exit 1
fi

echo "Input VCF: $VCF_FILE"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""

# Convert VCF to Relate format
echo "Running RelateFileFormats for chromosome X..."
RelateFileFormats \
    --mode ConvertFromVcf \
    --haps "${OUTPUT_PREFIX}.haps" \
    --sample "${OUTPUT_PREFIX}.sample" \
    -i "$VCF_FILE"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Conversion failed for chromosome X"
    exit 1
fi

echo ""
echo "✓ Conversion complete for chromosome X"
echo ""
echo "Output files:"
echo "  HAPS: ${OUTPUT_PREFIX}.haps"
echo "  SAMPLE: ${OUTPUT_PREFIX}.sample"
echo ""

# Show file sizes
ls -lh "${OUTPUT_PREFIX}.haps" "${OUTPUT_PREFIX}.sample"

echo ""
echo "=========================================="
echo "Important notes for chromosome X:"
echo "=========================================="
echo "- Males are hemizygous for most of chrX (except PAR regions)"
echo "- Relate handles this automatically if sex is specified in .sample file"
echo "- Use the 2/3 adjusted recombination map for chrX"
echo ""
echo "Next step: Run Relate inference with adjusted chrX genetic map"
