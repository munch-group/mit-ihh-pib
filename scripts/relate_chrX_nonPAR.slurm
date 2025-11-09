#!/bin/bash
#SBATCH --job-name=relate_chrX_nonPAR
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/relate_chrX_nonPAR_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/relate_chrX_nonPAR_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=normal

# Set the number of threads to match allocated CPUs
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Navigate to your Relate directory
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/relate/relate

echo "=========================================="
echo "Running Relate on chrX non-PAR region"
echo "Region: chrX:2,781,480-155,701,382 (non-PAR)"
echo "This excludes PAR1 and PAR2 to avoid variable haplotype counts"
echo "=========================================="
echo ""

# Run Relate using the parallel script for faster execution
# Using biallelic-only files to avoid duplicate position errors
./scripts/RelateParallel/RelateParallel.sh \
    --haps /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.haps \
    --sample /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.sample \
    --map /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR.relate.map \
    -m 0.83e-8 \
    -N 20000 \
    --seed 1 \
    --threads $SLURM_CPUS_PER_TASK \
    -o chrX_nonPAR_biallelic_output

echo ""
echo "Job completed at $(date)"
echo ""
echo "Output files:"
echo "  Trees: chrX_nonPAR_biallelic_output.anc"
echo "  Mutations: chrX_nonPAR_biallelic_output.mut"
echo ""
echo "Note: This analysis uses biallelic SNPs only (2,276,950 SNPs)"
echo "Covers non-PAR region only: chrX:2,781,480-155,701,382 (98% of chrX)"
echo "PAR regions excluded: PAR1 (10,001-2,781,479) and PAR2 (155,701,383-156,030,895)"
