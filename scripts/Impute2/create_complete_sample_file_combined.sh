#!/bin/bash
# Create complete sample file for all 3202 samples
# - Use known sex from panel for 2504 samples
# - Infer sex from heterozygosity for 698 related samples

set -euo pipefail

cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

VCF="/home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/chrX/chrX.vcf.gz"
PANEL="/home/vanbruggenmit/mit-ihh-pib/data/grch38/integrated_call_samples_v3.20130502.ALL.panel"
OUTPUT="/home/vanbruggenmit/mit-ihh-pib/data/grch38/raw/chrX/chrX.complete.sample"

echo "==================================================================="
echo "Creating complete sample file for 3202 samples"
echo "==================================================================="
echo ""

# Step 1: Get all samples from VCF
echo "Step 1: Extracting sample list from VCF..."
pixi run bcftools query -l "$VCF" > /tmp/all_samples.txt
total_samples=$(wc -l < /tmp/all_samples.txt)
echo "  Total samples in VCF: $total_samples"

# Step 2: Create associative array of known sexes from panel (2504 samples)
echo ""
echo "Step 2: Loading known sexes from panel file..."
declare -A known_sex
while IFS=$'\t' read -r sample pop super_pop sex; do
    if [[ "$sex" == "male" ]]; then
        known_sex[$sample]=1
    elif [[ "$sex" == "female" ]]; then
        known_sex[$sample]=2
    fi
done < <(tail -n +2 "$PANEL")

echo "  Known sexes from panel: ${#known_sex[@]} samples"

# Step 3: Identify samples that need sex inference (698 related samples)
echo ""
echo "Step 3: Identifying samples needing sex inference..."
unknown_samples=()
while read -r sample; do
    if [[ ! -v known_sex[$sample] ]]; then
        unknown_samples+=("$sample")
    fi
done < /tmp/all_samples.txt

echo "  Samples with known sex: ${#known_sex[@]}"
echo "  Samples needing inference: ${#unknown_samples[@]}"

# Step 4: Calculate heterozygosity for unknown samples only (FAST!)
if [[ ${#unknown_samples[@]} -gt 0 ]]; then
    echo ""
    echo "Step 4: Calculating heterozygosity for ${#unknown_samples[@]} samples..."
    
    # Create sample list for unknown samples
    printf '%s\n' "${unknown_samples[@]}" > /tmp/unknown_samples.txt
    
    # Calculate heterozygosity only for these samples
    pixi run bcftools view -S /tmp/unknown_samples.txt -r X:10000000-100000000 "$VCF" | \
    pixi run bcftools query -f '[%SAMPLE\t%GT\n]' | \
    awk '{
        sample = $1
        gt = $2
        total[sample]++
        if (gt == "0|1" || gt == "1|0") {
            het[sample]++
        }
    }
    END {
        for (s in total) {
            het_pct = (het[s] / total[s]) * 100
            # Threshold: <5% = male (1), ≥5% = female (2)
            sex = (het_pct < 5.0) ? 1 : 2
            print s, sex, het_pct
        }
    }' > /tmp/inferred_sex.txt
    
    echo "  Heterozygosity calculation complete"
    
    # Load inferred sexes
    declare -A inferred_sex
    declare -A het_values
    while read -r sample sex het_pct; do
        inferred_sex[$sample]=$sex
        het_values[$sample]=$het_pct
    done < /tmp/inferred_sex.txt
    
    # Count inferred males and females
    inferred_males=$(awk '$2 == 1' /tmp/inferred_sex.txt | wc -l)
    inferred_females=$(awk '$2 == 2' /tmp/inferred_sex.txt | wc -l)
    echo "  Inferred: $inferred_males males, $inferred_females females"
    
    # Show some examples
    echo ""
    echo "  Examples of inferred sex (first 10):"
    echo "  Sample       Sex    Het%"
    head -10 /tmp/inferred_sex.txt | awk '{printf "  %-12s %-6s %.2f%%\n", $1, ($2==1?"male":"female"), $3}'
fi

# Step 5: Create the sample file
echo ""
echo "Step 5: Writing sample file..."
echo "ID_1 ID_2 missing sex" > "$OUTPUT"
echo "0 0 0 D" >> "$OUTPUT"

males=0
females=0
from_panel=0
from_inference=0

while read -r sample; do
    # Use known sex if available, otherwise use inferred
    if [[ -v known_sex[$sample] ]]; then
        sex=${known_sex[$sample]}
        ((from_panel++))
    elif [[ -v inferred_sex[$sample] ]]; then
        sex=${inferred_sex[$sample]}
        ((from_inference++))
    else
        # Should not happen, but default to female if it does
        sex=2
        echo "  WARNING: No sex info for $sample, defaulting to female"
    fi
    
    # Count males and females
    if [[ $sex -eq 1 ]]; then
        ((males++))
    else
        ((females++))
    fi
    
    # Write to output
    echo "$sample $sample 0.0 $sex" >> "$OUTPUT"
done < /tmp/all_samples.txt

echo "  Complete!"

echo ""
echo "==================================================================="
echo "SAMPLE FILE CREATED"
echo "==================================================================="
echo "Output: $OUTPUT"
echo ""
echo "Sample counts:"
echo "  Total samples: $total_samples"
echo "  Males: $males ($from_panel from panel, $((males - from_panel)) inferred)"
echo "  Females: $females ($((from_panel - (males - from_inference))) from panel, $from_inference inferred)"
echo ""
echo "Sex source:"
echo "  From panel file: $from_panel samples"
echo "  Inferred from heterozygosity: $from_inference samples"
echo ""
echo "Expected haplotype counts:"
echo "  PAR regions: $((2 * total_samples)) (all diploid)"
echo "  Non-PAR (compact format): $((males + 2 * females))"
echo "  Non-PAR (with dummies): $((2 * total_samples))"
echo ""
echo "First 12 lines of sample file:"
head -12 "$OUTPUT"
echo ""
echo "Verifying sample file has correct number of lines:"
wc -l "$OUTPUT"
echo "  (Should be $((total_samples + 2)): 2 header + $total_samples samples)"