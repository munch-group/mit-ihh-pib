#!/bin/bash
#SBATCH --job-name=filter_chrX_biallelic
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/filter_chrX_biallelic_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/filter_chrX_biallelic_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=normal

echo "==========================================  "
echo "Filtering chrX VCF to biallelic SNPs only"
echo "Region: chrX:2,781,480-155,701,382 (non-PAR)"
echo "=========================================="
echo ""

# Activate pixi environment
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib
eval "$(pixi shell-hook)"

# Filter to region and SNPs only, split multiallelics, and remove duplicate positions
# -m -snps splits multiallelic sites into biallelic records
# -d snps removes duplicate SNP records (keeps first occurrence)
bcftools view -r chrX:2781480-155701382 -v snps \
  /faststorage/project/mit-ihh-pib/data/grch38/raw/chrX/chrX.vcf.gz \
  -O v | \
bcftools norm -m -snps -O v | \
bcftools norm -d snps -O z -o /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.vcf.gz

echo ""
echo "Indexing filtered VCF..."
bcftools index /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.vcf.gz

echo ""
echo "Filtering completed!"
echo "Output: /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.vcf.gz"

# Get SNP count
SNP_COUNT=$(bcftools view -H /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.vcf.gz | wc -l)
echo "Total biallelic SNPs: $SNP_COUNT"

echo ""
echo "Converting to Relate format..."
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/relate/relate

./bin/RelateFileFormats \
  --mode ConvertFromVcf \
  --haps /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.haps \
  --sample /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.sample \
  -i /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic

echo ""
echo "Verifying no duplicate positions..."
DUPLICATES=$(awk '{print $3}' /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.haps | sort -n | uniq -d | wc -l)
echo "Number of duplicate positions: $DUPLICATES"

if [ $DUPLICATES -eq 0 ]; then
    echo "✓ SUCCESS: No duplicate positions found!"
else
    echo "✗ WARNING: Still found $DUPLICATES duplicate positions"
fi

echo ""
echo "Job completed at $(date)"
