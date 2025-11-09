#!/bin/bash
#SBATCH --job-name=check_vcf
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/check_vcf_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/check_vcf_%j.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=normal

cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib
eval "$(pixi shell-hook)"

echo "Checking VCF at position 2783144:"
bcftools view -H /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.vcf.gz | awk '$2 == 2783144' | cut -f1-5

echo ""
echo "Checking how many times position 2783144 appears:"
bcftools view -H /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.vcf.gz | awk '$2 == 2783144' | wc -l

echo ""
echo "Checking a few duplicate positions from VCF:"
bcftools view -H /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.vcf.gz | awk '{print $2}' | sort -n | uniq -d | head -5 | while read pos; do
    echo "Position $pos:"
    bcftools view -H /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.vcf.gz | awk -v p=$pos '$2 == p' | cut -f1-5
    echo ""
done

echo "Total SNPs in VCF:"
bcftools view -H /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.vcf.gz | wc -l

echo ""
echo "Checking for multiallelic sites (should be 0):"
bcftools view -H /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.vcf.gz | awk '{split($5,a,","); if(length(a)>1) print $1,$2,$4,$5}' | wc -l
