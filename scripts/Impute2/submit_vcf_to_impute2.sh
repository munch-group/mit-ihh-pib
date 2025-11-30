#!/bin/bash
#SBATCH --job-name=vcf2impute2
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/vcf2impute2_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/vcf2impute2_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal

# Run the conversion script with pixi environment
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib

# Execute with pixi
pixi run bash scripts/Impute2/run_vcf_to_impute2_conversion.sh
