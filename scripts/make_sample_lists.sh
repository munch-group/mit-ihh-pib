#!/usr/bin/env bash
set -euo pipefail
DATA_ROOT=${1:-/home/vanbruggenmit/mit-ihh-pib/data}

mkdir -p "${DATA_ROOT}/panels"

wget -c -P "${DATA_ROOT}/panels" \
  https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# Example populations (edit or duplicate these lines as needed)
awk '$2=="GBR"{print $1}' "${DATA_ROOT}/panels/integrated_call_samples_v3.20130502.ALL.panel" > "${DATA_ROOT}/panels/GBR.samples"
awk '$2=="YRI"{print $1}' "${DATA_ROOT}/panels/integrated_call_samples_v3.20130502.ALL.panel" > "${DATA_ROOT}/panels/YRI.samples"
