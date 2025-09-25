#!/usr/bin/env bash
set -euo pipefail
DATA_ROOT=${1:-/home/vanbruggenmit/mit-ihh-pib/data}

# Space-separated population codes; edit as needed.
POPS=${POPS:-"GBR YRI"}

# change CHRS from "22" to "{1..22}" after testing
CHRS=${CHRS:-22}

mkdir -p "${DATA_ROOT}/subsets"

for chr in ${CHRS}; do
  vcf="${DATA_ROOT}/vcfs/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  for POP in ${POPS}; do
    pixi run bcftools view -S "${DATA_ROOT}/panels/${POP}.samples" -m2 -M2 -v snps \
      -Oz -o "${DATA_ROOT}/subsets/${POP}.chr${chr}.snps.vcf.gz" "${vcf}"
    pixi run tabix -p vcf "${DATA_ROOT}/subsets/${POP}.chr${chr}.snps.vcf.gz"
  done
done
