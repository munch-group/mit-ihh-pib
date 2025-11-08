#!/usr/bin/env bash
#SBATCH --job-name=vcf_to_relate
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/vcf_to_relate_%A_%a.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/vcf_to_relate_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22
#SBATCH --partition=normal

# Convert VCF files to Relate format for selection analysis verification
# This script handles autosomes (chromosomes 1-22)

# Set up Relate path
RELATE_BIN="/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/software/relate/bin"
export PATH="$RELATE_BIN:$PATH"

# Check if Relate is available
if ! command -v RelateFileFormats &> /dev/null; then
    echo "ERROR: RelateFileFormats not found"
    echo "Please run: sbatch scripts/install_relate.sh"
    exit 1
fi

# Get chromosome from SLURM array task ID
CHR=${SLURM_ARRAY_TASK_ID}

echo "=========================================="
echo "Converting chromosome ${CHR} to Relate format"
echo "=========================================="

# Paths
VCF_DIR="/faststorage/project/mit-ihh-pib/data/grch38/raw/chr${CHR}"
VCF_FILE="$VCF_DIR/chr${CHR}.vcf.gz"
OUTPUT_DIR="/faststorage/project/mit-ihh-pib/data/grch38/relate"
OUTPUT_PREFIX="$OUTPUT_DIR/chr${CHR}"

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
echo "Running RelateFileFormats..."
RelateFileFormats \
    --mode ConvertFromVcf \
    --haps "${OUTPUT_PREFIX}.haps" \
    --sample "${OUTPUT_PREFIX}.sample" \
    -i "$VCF_FILE"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Conversion failed for chromosome ${CHR}"
    exit 1
fi

echo ""
echo "✓ Conversion complete for chromosome ${CHR}"
echo ""
echo "Output files:"
echo "  HAPS: ${OUTPUT_PREFIX}.haps"
echo "  SAMPLE: ${OUTPUT_PREFIX}.sample"
echo ""

# Show file sizes
ls -lh "${OUTPUT_PREFIX}.haps" "${OUTPUT_PREFIX}.sample"

echo ""
echo "Next step: Prepare genetic map and run Relate inference"
