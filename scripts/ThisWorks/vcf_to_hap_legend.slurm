#!/usr/bin/env bash
#SBATCH --job-name=vcf2hap_chr%a
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/vcf2hap_chr%a_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/vcf2hap_chr%a_%j.err
#SBATCH --array=1-20
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=normal

# Get chromosome number from array task ID
CHR=$SLURM_ARRAY_TASK_ID

DATA_ROOT="/home/vanbruggenmit/mit-ihh-pib/data/grch38"
RAW_DIR="$DATA_ROOT/raw/chr${CHR}"
WORK_DIR="$DATA_ROOT/work/chr${CHR}"

mkdir -p "$RAW_DIR" "$WORK_DIR"

echo "=========================================="
echo "Downloading and converting chromosome $CHR"
echo "=========================================="

# VCF URL from config.yaml
VCF_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz"
VCF_FILE="$RAW_DIR/chr${CHR}.vcf.gz"
VCF_INDEX="$VCF_FILE.tbi"

# Download VCF if not exists
if [[ ! -f "$VCF_FILE" ]]; then
    echo "Downloading VCF for chr${CHR}..."
    wget -O "$VCF_FILE" "$VCF_URL"
    
    if [[ $? -ne 0 ]]; then
        echo "ERROR: Failed to download VCF"
        exit 1
    fi
    
    # Download index
    echo "Downloading VCF index..."
    wget -O "$VCF_INDEX" "${VCF_URL}.tbi"
else
    echo "VCF already exists: $VCF_FILE"
fi

# Convert VCF to HAP/LEGEND using bcftools
echo ""
echo "Converting VCF to HAP/LEGEND format..."

# Convert to HAP (haplotype matrix)
bcftools query -f '[%GT\t]\n' "$VCF_FILE" | \
    sed 's/|/ /g' | sed 's/\t/ /g' | \
    gzip > "$WORK_DIR/chr${CHR}.impute.hap.gz"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to create HAP file"
    exit 1
fi

# Convert to LEGEND (variant information)
echo "id position a0 a1" > "$WORK_DIR/chr${CHR}.impute.legend"
bcftools query -f '%ID\t%POS\t%REF\t%ALT\n' "$VCF_FILE" >> "$WORK_DIR/chr${CHR}.impute.legend"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to create LEGEND file"
    exit 1
fi

gzip "$WORK_DIR/chr${CHR}.impute.legend"

echo ""
echo "✓ Conversion complete for chr${CHR}"
echo "  HAP: $WORK_DIR/chr${CHR}.impute.hap.gz"
echo "  LEGEND: $WORK_DIR/chr${CHR}.impute.legend.gz"

# Show file sizes
ls -lh "$WORK_DIR/chr${CHR}.impute."*