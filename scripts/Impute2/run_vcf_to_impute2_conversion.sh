#!/bin/bash
set -euo pipefail

cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

WORK_DIR="/home/vanbruggenmit/mit-ihh-pib/data/grch38/work/chrX_impute2_conversion"
VCF="/home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/X.vcf.gz"
SAMPLE_FILE="/home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/chrX/chrX.complete.sample"
PED_FILE="/home/vanbruggenmit/mit-ihh-pib/data/grch38/integrated_call_samples_v3.20200731.ALL.ped"

cd "$WORK_DIR"

echo "==================================================================="
echo "Step 1: Creating sample table in format the Perl script expects"
echo "==================================================================="

# The Perl script expects: sample_id  pop_id  group_id  sex
# We have PED file with: FamilyID IndividualID PaternalID MaternalID Gender Phenotype Population ...
# We need to create: sample_id pop_id group_id sex

tail -n +2 "$PED_FILE" | \
awk 'BEGIN {OFS="\t"} {
    sample_id = $2
    pop_id = $7        # Population column
    group_id = $7      # Use population as group for now
    sex = $5           # Gender: 1=male, 2=female
    print sample_id, pop_id, group_id, sex
}' > chrX.sample_table.txt

echo "Sample table created with $(wc -l < chrX.sample_table.txt) samples"
echo ""
echo "First 10 lines:"
head -10 chrX.sample_table.txt

echo ""
echo "==================================================================="
echo "Step 2: Running IMPUTE2 conversion (FULL chrX including PAR)"
echo "==================================================================="
echo ""
echo "This will create:"
echo "  - chrX.impute.legend.gz"
echo "  - chrX.impute.hap.gz"
echo "  - chrX.impute.sample_list"
echo ""

# Run the Perl script for FULL chrX (PAR + non-PAR)
perl /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/scripts/Impute2/vcf2impute_legend_haps.pl \
    -vcf "$VCF" \
    -leghap chrX.impute \
    -chr X \
    -samp_tab chrX.sample_table.txt

echo ""
echo "==================================================================="
echo "Step 3: Verifying output files"
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