#!/usr/bin/bash
#SBATCH --job-name=count_hets
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/count_hets_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/count_hets_%j.err

cd /home/vanbruggenmit/mit-ihh-pib/data/grch38/work/chrX_impute2_conversion

echo "Scanning chrX.nonPAR.vcf.gz for positions where males have heterozygous genotypes..."
echo ""

# Scan the entire nonPAR VCF and count variants with male hets
zcat chrX.nonPAR.vcf.gz | grep -v '^#' | awk '
BEGIN {
  problematic_count = 0
  total_count = 0
}
{
  total_count++
  het_count = 0

  # Check all sample columns (starting from column 10)
  for(i=10; i<=NF; i++) {
    if($i == "0|1" || $i == "1|0") {
      het_count++
    }
  }

  if(het_count > 0) {
    problematic_count++
    # Print first 50 problematic positions
    if(problematic_count <= 50) {
      print "Position", $2, "has", het_count, "male(s) with het genotypes"
    }
  }

  # Progress indicator every 100k variants
  if(total_count % 100000 == 0) {
    print "Processed", total_count, "variants, found", problematic_count, "problematic so far..." > "/dev/stderr"
  }
}
END {
  print ""
  print "================================"
  print "SUMMARY:"
  print "Total variants scanned:", total_count
  print "Variants with male hets:", problematic_count
  print "Percentage:", (problematic_count/total_count*100), "%"
  print "================================"
}
'
