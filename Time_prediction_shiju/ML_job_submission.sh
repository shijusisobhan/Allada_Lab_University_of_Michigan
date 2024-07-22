#!/bin/bash
#SBATCH --job-name=ML
#SBATCH --nodes=2 
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=128g
#SBATCH --time=14:15:00
#SBATCH --account=rallada0
#SBATCH --partition=largemem
#SBATCH --mail-user=shijusis@umich.edu
#SBATCH --mail-type=END 
my_job_header
module load R
echo "Running from $(pwd)"

R CMD BATCH --no-restore --no-save --quiet Gene_ratio_ML_Bash.R ML_run_report.out