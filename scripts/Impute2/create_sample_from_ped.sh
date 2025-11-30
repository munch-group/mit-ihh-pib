#!/bin/bash
# Create complete sample file using PED file sex information

set -euo pipefail

cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

VCF="/home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/chrX/chrX.vcf.gz"
PED="/home/vanbruggenmit/mit-ihh-pib/data/grch38/integrated_call_samples_v3.20200731.ALL.ped"
OUTPUT="/home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/chrX/chrX.complete.sample"

echo "==================================================================="
echo "Creating complete sample file from PED file"
echo "==================================================================="
echo ""

# Step 1: Get all samples from VCF
echo "Step 1: Extracting sample list from VCF..."
pixi run bcftools query -l "$VCF" > /tmp/all_samples.txt
total_samples=$(wc -l < /tmp/all_samples.txt)
echo "  Total samples in VCF: $total_samples"

# Step 2: Load sex information from PED file
echo ""
echo "Step 2: Loading sex information from PED file..."
declare -A sample_sex

while IFS=$'\t' read -r fam_id ind_id pat_id mat_id gender pheno rest; do
    # Skip header line
    if [[ "$ind_id" == "Individual ID" ]]; then
        continue
    fi
    
    # Store sex (column 5: Gender, 1=male, 2=female)
    if [[ "$gender" == "1" ]] || [[ "$gender" == "2" ]]; then
        sample_sex[$ind_id]=$gender
    fi
done < "$PED"

echo "  Loaded sex info for ${#sample_sex[@]} samples from PED file"

# Step 3: Create sample file
echo ""
echo "Step 3: Creating IMPUTE2 sample file..."
echo "ID_1 ID_2 missing sex" > "$OUTPUT"
echo "0 0 0 D" >> "$OUTPUT"

males=0
females=0
found=0
missing=0
missing_samples=()

while read -r sample; do
    if [[ -v sample_sex[$sample] ]]; then
        sex=${sample_sex[$sample]}
        ((found++))
    else
        # Default to female (2) if not found - safer for diploid assumption
        sex=2
        ((missing++))
        missing_samples+=("$sample")
    fi
    
    # Count
    if [[ $sex -eq 1 ]]; then
        ((males++))
    else
        ((females++))
    fi
    
    # Write
    echo "$sample $sample 0.0 $sex" >> "$OUTPUT"
done < /tmp/all_samples.txt

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
echo "  Missing from PED: $missing"

if [[ $missing -gt 0 ]]; then
    echo ""
    echo "Missing samples (defaulted to female):"
    printf '  %s\n' "${missing_samples[@]}" | head -20
    if [[ $missing -gt 20 ]]; then
        echo "  ... and $((missing - 20)) more"
    fi
fi

echo ""
echo "Expected haplotype counts:"
echo "  PAR regions: $((2 * total_samples))"
echo "  Non-PAR (compact format): $((males + 2 * females))"
echo "  Non-PAR (with dummies): $((2 * total_samples))"
echo ""
echo "First 12 lines of sample file:"
head -12 "$OUTPUT"
echo ""
echo "Verification:"
wc -l "$OUTPUT"
echo "  (Should be $((total_samples + 2)): 2 header + $total_samples samples)"