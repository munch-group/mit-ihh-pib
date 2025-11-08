#!/bin/bash
#SBATCH --job-name=relate_chrX
#SBATCH --output=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/relate_chrX_%j.out
#SBATCH --error=/home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/logs/relate_chrX_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=normal

# Set the number of threads to match allocated CPUs
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Navigate to your Relate directory
cd /home/vanbruggenmit/mit-ihh-pib/people/vanbruggenmit/mit-ihh-pib/relate/relate

# Run Relate using the parallel script for faster execution
./scripts/RelateParallel/RelateParallel.sh \
    --haps /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX/chrX.haps \
    --sample /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX/chrX.sample \
    --map /home/vanbruggenmit/mit-ihh-pib/data/relate/chrX/chrX.relate.map \
    -m 0.83e-8 \
    -N 20000 \
    --seed 1 \
    --threads $SLURM_CPUS_PER_TASK \
    -o chrX_output

echo "Job completed at $(date)"
