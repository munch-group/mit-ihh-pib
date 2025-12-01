#!/bin/bash
set -euo pipefail

cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

WORK_DIR="/home/vanbruggenmit/mit-ihh-pib/data/grch38/work/chrX_impute2_conversion"
VCF="/home/vanbruggenmit/mit-ihh-pib/data/grch38/work/chrX_impute2_conversion/X.normalized.vcf.gz"
SAMPLE_FILE="/home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/chrX/chrX.complete.sample"
PED_FILE="/home/vanbruggenmit/mit-ihh-pib/data/grch38/integrated_call_samples_v3.20200731.ALL.ped"

cd "$WORK_DIR"

echo "==================================================================="
echo "Step 1: Creating sample table in VCF order"
echo "==================================================================="

# CRITICAL: Sample table must be in the same order as samples in the VCF
# The Perl script expects: sample_id  pop_id  group_id  sex

# First, get the sample list from VCF (in VCF order)
bcftools query -l "$VCF" > chrX.vcf_samples.txt

# Then, create sample table by looking up each VCF sample in the PED file
while read sample_id; do
    # Look up this sample in the PED file
    tail -n +2 "$PED_FILE" | awk -v sample="$sample_id" '
        $2 == sample {
            pop_id = $7        # Population column
            group_id = $7      # Use population as group
            sex = $5           # Gender: 1=male, 2=female
            print sample, pop_id, group_id, sex
        }
    '
done < chrX.vcf_samples.txt > chrX.sample_table.txt

echo "Sample table created with $(wc -l < chrX.sample_table.txt) samples (matching VCF order)"
echo ""
echo "First 10 lines:"
head -10 chrX.sample_table.txt

echo ""
echo "==================================================================="
echo "Step 2: Indexing VCF file if needed"
echo "==================================================================="
echo ""

# Check if index exists, create if needed
if [[ ! -f "${VCF}.tbi" ]]; then
    echo "Index not found. Creating tabix index for VCF file..."
    tabix -p vcf "$VCF"
    echo "✓ Index created: ${VCF}.tbi"
else
    echo "✓ Index already exists: ${VCF}.tbi"
fi

echo ""
echo "==================================================================="
echo "Step 3: Splitting chrX VCF by region (PAR1, non-PAR, PAR2)"
echo "==================================================================="
echo ""

# GRCh38 coordinates for chrX regions
PAR1_START=60001
PAR1_END=2781479
NONPAR_START=2781480  # Start immediately after PAR1
NONPAR_END=155701383  # End right before PAR2
PAR2_START=155701384
PAR2_END=156030895

echo "Creating region-specific VCF files..."
echo "  PAR1:    X:${PAR1_START}-${PAR1_END}"
echo "  non-PAR: X:${NONPAR_START}-${NONPAR_END}"
echo "  PAR2:    X:${PAR2_START}-${PAR2_END}"
echo ""

# Extract PAR1
bcftools view -r X:${PAR1_START}-${PAR1_END} "$VCF" -Oz -o chrX.PAR1.vcf.gz
tabix -p vcf chrX.PAR1.vcf.gz

# Extract non-PAR region
echo "  Extracting non-PAR region..."
bcftools view -r X:${NONPAR_START}-${NONPAR_END} "$VCF" -Oz -o chrX.nonPAR.vcf.gz
tabix -p vcf chrX.nonPAR.vcf.gz

echo "  Note: Using normalized VCF where male genotypes have been converted:"
echo "    - Homozygous diploid (0|0, 1|1) → Haploid (0, 1)"
echo "    - Heterozygous diploid (0|1, 1|0) → Missing (.)"

# Extract PAR2
bcftools view -r X:${PAR2_START}-${PAR2_END} "$VCF" -Oz -o chrX.PAR2.vcf.gz
tabix -p vcf chrX.PAR2.vcf.gz

echo "VCF files created:"
ls -lh chrX.PAR*.vcf.gz chrX.nonPAR.vcf.gz

echo ""
echo "==================================================================="
echo "Step 4: Running IMPUTE2 conversion for each region"
echo "==================================================================="
echo ""

# Convert PAR1 (diploid in both males and females)
echo "Converting PAR1..."
perl /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/Impute2/vcf2impute_legend_haps.pl \
    -vcf chrX.PAR1.vcf.gz \
    -leghap chrX.PAR1.impute \
    -chr X \
    -samp_tab chrX.sample_table.txt

# Convert non-PAR (hemizygous in males)
echo ""
echo "Converting non-PAR (with -chrX_nonpar flag)..."
perl /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/Impute2/vcf2impute_legend_haps.pl \
    -vcf chrX.nonPAR.vcf.gz \
    -leghap chrX.nonPAR.impute \
    -chr X \
    -chrX_nonpar \
    -samp_tab chrX.sample_table.txt

# Convert PAR2 (diploid in both males and females)
echo ""
echo "Converting PAR2..."
perl /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/Impute2/vcf2impute_legend_haps.pl \
    -vcf chrX.PAR2.vcf.gz \
    -leghap chrX.PAR2.impute \
    -chr X \
    -samp_tab chrX.sample_table.txt

echo ""
echo "==================================================================="
echo "Step 5: Combining regions into single files"
echo "==================================================================="
echo ""

# Combine legend files (skip header from non-PAR and PAR2)
echo "Combining legend files..."
zcat chrX.PAR1.impute.legend.gz > chrX.impute.legend
zcat chrX.nonPAR.impute.legend.gz | tail -n +2 >> chrX.impute.legend
zcat chrX.PAR2.impute.legend.gz | tail -n +2 >> chrX.impute.legend
gzip chrX.impute.legend

# Combine haplotype files
echo "Combining haplotype files..."
zcat chrX.PAR1.impute.hap.gz > chrX.impute.hap
zcat chrX.nonPAR.impute.hap.gz >> chrX.impute.hap
zcat chrX.PAR2.impute.hap.gz >> chrX.impute.hap
gzip chrX.impute.hap

# Use sample list from PAR1 (should be identical across all regions)
cp chrX.PAR1.impute.sample_list chrX.impute.sample_list

echo "Combined files created:"
ls -lh chrX.impute.legend.gz chrX.impute.hap.gz chrX.impute.sample_list

echo ""
echo "==================================================================="
echo "Step 6: Verifying output files"
echo "==================================================================="

if [[ -f chrX.impute.legend.gz ]] && [[ -f chrX.impute.hap.gz ]]; then
    echo "✓ Files created successfully!"
    echo ""
    ls -lh chrX.impute.*
    
    echo ""
    echo "=== Checking haplotype file ==="
    echo "First variant (should be in PAR1):"
    zcat chrX.impute.hap.gz | head -1 | awk '{print "Fields:", NF}'
    
    echo ""
    echo "Checking consistency (first 1000 lines):"
    zcat chrX.impute.hap.gz | head -1000 | awk '{print NF}' | sort | uniq -c
    
    echo ""
    echo "=== Checking legend file ==="
    echo "Header:"
    zcat chrX.impute.legend.gz | head -1
    echo ""
    echo "First 3 variants:"
    zcat chrX.impute.legend.gz | head -4 | tail -3
    
else
    echo "✗ ERROR: Output files not created!"
    exit 1
fi

echo ""
echo "==================================================================="
echo "COMPLETE!"
echo "==================================================================="
echo "Output location: $WORK_DIR"
echo ""
echo "Next steps:"
echo "1. Check that field counts are consistent"
echo "2. Convert to hapbin format"
echo "3. Run iHS"