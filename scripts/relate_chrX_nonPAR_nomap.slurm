#!/bin/bash
#SBATCH --job-name=relate_chrX_nomap
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/relate_chrX_nomap_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/relate_chrX_nomap_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=normal

# Set the number of threads to match allocated CPUs
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Navigate to your Relate directory
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/relate/relate

echo "=========================================="
echo "Running Relate WITHOUT genetic map on chrX non-PAR region"
echo "Region: chrX:3,533,229-154,781,072"
echo "Will use constant recombination rate estimated from mutation rate"
echo "=========================================="
echo ""

# Try running Relate without genetic map
# Relate should be able to work with just distances and mutation rate
./bin/Relate \
    --mode All \
    --haps /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.haps \
    --sample /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.sample \
    --dist /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.dist \
    -m 0.83e-8 \
    -N 20000 \
    --seed 1 \
    -o chrX_nonPAR_nomap_output

echo ""
echo "Job completed at $(date)"
