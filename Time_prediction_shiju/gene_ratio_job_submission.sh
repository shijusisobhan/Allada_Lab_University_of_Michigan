#!/bin/bash
#SBATCH --job-name=Gene_ratio
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=128g
#SBATCH --time=04:15:00
#SBATCH --account=rallada0
#SBATCH --partition=standard
#SBATCH --mail-user=shijusis@umich.edu
#SBATCH --mail-type=END 
my_job_header
module load R
echo "Running from $(pwd)"

R CMD BATCH --no-restore --no-save --quiet Gene-ratio_function_bash.R MAB_z_score.out