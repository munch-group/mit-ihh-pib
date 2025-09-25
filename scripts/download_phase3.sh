#!/usr/bin/env bash
set -euo pipefail
DATA_ROOT=${1:-/home/vanbruggenmit/mit-ihh-pib/data}

prefix="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr"
suffix=".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

mkdir -p "${DATA_ROOT}/vcfs"

# change CHRS from "22" to "{1..22}" after testing
CHRS=${CHRS:-22}

for chr in ${CHRS}; do
  wget -c -P "${DATA_ROOT}/vcfs" "${prefix}${chr}${suffix}"
  wget -c -P "${DATA_ROOT}/vcfs" "${prefix}${chr}${suffix}.tbi"
done
