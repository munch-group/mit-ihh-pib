#!/bin/bash
# Create complete sample file using PED file - FAST VERSION

set -euo pipefail

cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

VCF="/home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/chrX/chrX.vcf.gz"
PED="/home/vanbruggenmit/mit-ihh-pib/data/grch38/integrated_call_samples_v3.20200731.ALL.ped"
OUTPUT="/home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/chrX/chrX.complete.sample"

echo "==================================================================="
echo "Creating complete sample file from PED file (FAST)"
echo "==================================================================="
echo ""

# Step 1: Get all samples from VCF
echo "Step 1: Extracting sample list from VCF..."
pixi run bcftools query -l "$VCF" > /tmp/vcf_samples.txt
total_samples=$(wc -l < /tmp/vcf_samples.txt)
echo "  Total samples in VCF: $total_samples"

# Step 2: Extract sample ID and sex from PED file
echo ""
echo "Step 2: Extracting sex from PED file..."
tail -n +2 "$PED" | awk '{print $2, $5}' | sort -k1,1 > /tmp/ped_sex.txt
echo "  Extracted sex for $(wc -l < /tmp/ped_sex.txt) samples"

# Step 3: Sort VCF samples
echo ""
echo "Step 3: Preparing sample list..."
sort /tmp/vcf_samples.txt > /tmp/vcf_samples_sorted.txt

# Step 4: Join VCF samples with PED sex info
echo ""
echo "Step 4: Joining sample and sex information..."
join -a 1 /tmp/vcf_samples_sorted.txt /tmp/ped_sex.txt | \
    awk '{
        sample = $1
        sex = ($2 == "") ? 2 : $2  # Default to female (2) if missing
        print sample, sample, "0.0", sex
    }' > /tmp/samples_with_sex.txt

# Step 5: Count males and females
males=$(awk '$4 == 1' /tmp/samples_with_sex.txt | wc -l)
females=$(awk '$4 == 2' /tmp/samples_with_sex.txt | wc -l)
found=$(join /tmp/vcf_samples_sorted.txt /tmp/ped_sex.txt | wc -l)
missing=$((total_samples - found))

echo "  Males: $males"
echo "  Females: $females"
echo "  Found in PED: $found"
echo "  Missing from PED: $missing"

# Step 6: Create IMPUTE2 sample file
echo ""
echo "Step 5: Writing sample file..."
echo "ID_1 ID_2 missing sex" > "$OUTPUT"
echo "0 0 0 D" >> "$OUTPUT"
cat /tmp/samples_with_sex.txt >> "$OUTPUT"

echo ""
echo "==================================================================="
echo "COMPLETE!"
echo "==================================================================="
echo "Sample file: $OUTPUT"
echo ""
echo "Summary:"
echo "  Total samples: $total_samples"
echo "  Males: $males"
echo "  Females: $females"
echo "  Found in PED: $found"
echo "  Missing from PED: $missing (defaulted to female)"
echo ""
echo "Expected haplotype counts:"
echo "  PAR regions: $((2 * total_samples))"
echo "  Non-PAR (compact format): $((males + 2 * females))"
echo "  Non-PAR (with dummies): $((2 * total_samples))"
echo ""
echo "First 12 lines of sample file:"
head -12 "$OUTPUT"
echo ""
echo "Last 5 lines:"
tail -5 "$OUTPUT"
echo ""
echo "Verification:"
echo "Lines in file: $(wc -l < "$OUTPUT")"
echo "Expected: $((total_samples + 2))"