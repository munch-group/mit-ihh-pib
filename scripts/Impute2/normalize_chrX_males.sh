#!/usr/bin/bash
#SBATCH --job-name=normalize_chrX
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/normalize_chrX_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/normalize_chrX_%j.err

cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

echo "================================================================="
echo "Normalizing chromosome X VCF for Eagle2 phasing output"
echo "================================================================="
echo ""
echo "Converting male genotypes in non-PAR region (2,781,480-155,701,383):"
echo "  - Homozygous diploid (0|0, 1|1) → Haploid (0, 1)"
echo "  - Heterozygous diploid (0|1, 1|0) → Missing (.)"
echo "  - Female genotypes unchanged"
echo ""

VCF_IN="/home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/X.vcf.gz"
VCF_OUT="/home/vanbruggenmit/mit-ihh-pib/data/grch38/work/chrX_impute2_conversion/X.normalized.vcf.gz"
PED_FILE="/faststorage/project/mit-ihh-pib/data/grch38/integrated_call_samples_v3.20200731.ALL.ped"

# Extract male sample IDs that are actually in the VCF
echo "Step 1: Creating male samples file..."
# Get VCF samples
pixi run bcftools query -l "$VCF_IN" > /tmp/vcf_samples.txt
# Get male samples from PED
awk '$5 == 1 {print $2}' "$PED_FILE" > /tmp/ped_males.txt
# Intersect: only males that are in both PED and VCF
grep -Fx -f /tmp/vcf_samples.txt /tmp/ped_males.txt > /tmp/male_samples.txt
n_males=$(wc -l < /tmp/male_samples.txt)
echo "  Found $n_males male samples (present in both PED and VCF)"
rm -f /tmp/vcf_samples.txt /tmp/ped_males.txt
echo ""

# Process PAR1 (keep as-is, diploid for all)
echo "Step 2: Processing PAR1 region (60,001-2,781,479)..."
pixi run bcftools view -r X:60001-2781479 "$VCF_IN" -Oz -o /tmp/X.PAR1.vcf.gz
pixi run tabix -p vcf /tmp/X.PAR1.vcf.gz
echo "  ✓ PAR1 extracted (no changes needed)"
echo ""

# Process non-PAR with male genotype conversion
echo "Step 3: Processing non-PAR region (2,781,480-155,701,383)..."
echo "  Extracting region and processing males separately..."

# Step 3a: Extract males from non-PAR and convert genotypes
pixi run bcftools view -r X:2781480-155701383 -S /tmp/male_samples.txt "$VCF_IN" -Ou | \
pixi run bcftools +setGT -Ou -- -t q -i 'GT="0|0"' -n c:0 | \
pixi run bcftools +setGT -Ou -- -t q -i 'GT="1|1"' -n c:1 | \
pixi run bcftools +setGT -Oz -o /tmp/X.nonPAR.males.vcf.gz -- -t q -i 'GT="het"' -n c:.
pixi run tabix -p vcf /tmp/X.nonPAR.males.vcf.gz

echo "  Converting male genotypes (0|0→0, 1|1→1, het→.)..."

# Step 3b: Extract females from non-PAR (unchanged)
pixi run bcftools query -l "$VCF_IN" | grep -v -F -f /tmp/male_samples.txt > /tmp/female_samples.txt
pixi run bcftools view -r X:2781480-155701383 -S /tmp/female_samples.txt "$VCF_IN" -Oz -o /tmp/X.nonPAR.females.vcf.gz
pixi run tabix -p vcf /tmp/X.nonPAR.females.vcf.gz

echo "  Merging males and females back together..."

# Step 3c: Get original sample order
pixi run bcftools query -l "$VCF_IN" > /tmp/original_sample_order.txt

# Merge males and females
pixi run bcftools merge --force-samples /tmp/X.nonPAR.males.vcf.gz /tmp/X.nonPAR.females.vcf.gz -Ou | \
pixi run bcftools view -S /tmp/original_sample_order.txt -Oz -o /tmp/X.nonPAR.normalized.vcf.gz

pixi run tabix -p vcf /tmp/X.nonPAR.normalized.vcf.gz
echo "  ✓ Non-PAR normalized"
echo ""

# Process PAR2 (keep as-is, diploid for all)
echo "Step 4: Processing PAR2 region (155,701,384-156,030,895)..."
pixi run bcftools view -r X:155701384-156030895 "$VCF_IN" -Oz -o /tmp/X.PAR2.vcf.gz
pixi run tabix -p vcf /tmp/X.PAR2.vcf.gz
echo "  ✓ PAR2 extracted (no changes needed)"
echo ""

# Concatenate all regions
echo "Step 5: Concatenating PAR1, normalized non-PAR, and PAR2..."
pixi run bcftools concat \
  /tmp/X.PAR1.vcf.gz \
  /tmp/X.nonPAR.normalized.vcf.gz \
  /tmp/X.PAR2.vcf.gz \
  -Oz -o "$VCF_OUT"

pixi run tabix -p vcf "$VCF_OUT"
echo "  ✓ Normalized VCF created"
echo ""

# Cleanup
echo "Step 6: Cleaning up temporary files..."
rm -f /tmp/X.*.vcf.gz /tmp/X.*.vcf.gz.tbi /tmp/male_samples.txt /tmp/female_samples.txt
echo "  ✓ Cleanup complete"
echo ""

echo "================================================================="
echo "COMPLETED!"
echo "Normalized VCF: $VCF_OUT"
echo "================================================================="

# Verify the output
echo ""
echo "Verification: Checking genotype patterns at position 2,781,514..."
pixi run bcftools view -H -r X:2781514 "$VCF_OUT" | awk '
{
  hap_0=0; hap_1=0; dip_hom=0; dip_het=0; missing=0
  for(i=10;i<=NF;i++) {
    if($i==".") missing++
    else if($i=="0") hap_0++
    else if($i=="1") hap_1++
    else if($i=="0|0"||$i=="1|1") dip_hom++
    else if($i=="0|1"||$i=="1|0") dip_het++
  }
  print "After normalization:"
  print "  Haploid 0:", hap_0
  print "  Haploid 1:", hap_1
  print "  Diploid homozygous:", dip_hom
  print "  Diploid heterozygous:", dip_het
  print "  Missing:", missing
  print ""
  print "Expected: ~1600 haploid (males), ~1600 diploid (females), ~700 missing (male hets converted)"
}
'
