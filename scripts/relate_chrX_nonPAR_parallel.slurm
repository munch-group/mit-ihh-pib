#!/bin/bash
#SBATCH --job-name=relate_chrX_nonPAR_parallel
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/relate_chrX_nonPAR_parallel_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/relate_chrX_nonPAR_parallel_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=normal

# Set the number of threads to match allocated CPUs
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Navigate to your Relate directory
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/relate/relate

echo "=========================================="
echo "Running RelateParallel on chrX non-PAR region"
echo "Region: chrX:3,533,229-154,781,072"
echo "Using biallelic SNPs only (2,252,381 SNPs)"
echo "=========================================="
echo ""

# Run Relate using the proper parallel script
# Using sparse genetic map (23,976 positions) - Relate will interpolate
# Using dist file for precise SNP distances
./scripts/RelateParallel/RelateParallel.sh \
    --haps /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.haps \
    --sample /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.sample \
    --map /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR.relate.map \
    --dist /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX_nonPAR/chrX_nonPAR_biallelic.dist \
    -m 0.83e-8 \
    -N 20000 \
    --seed 1 \
    --threads $SLURM_CPUS_PER_TASK \
    -o chrX_nonPAR_parallel_output

echo ""
echo "Job completed at $(date)"
echo ""
echo "Output files:"
echo "  Trees: chrX_nonPAR_parallel_output.anc"
echo "  Mutations: chrX_nonPAR_parallel_output.mut"
