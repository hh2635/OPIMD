#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=4GB
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --account=torch_pr_226_chemistry

echo "Job started at $(date)"

module purge 
module load anaconda3/2025.06

make clean 
make 
make run
make clean

python ../analysis/Wigner_exact.py
python ../analysis/Wigner_marginal.py

echo "Job finished at $(date)"