#!/bin/bash
#SBATCH --job-name=FLIC-1
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=32g
#SBATCH --time=01:15:00
#SBATCH --account=rallada0
#SBATCH --partition=standard
#SBATCH --mail-user=shijusis@umich.edu
#SBATCH --mail-type=END 
my_job_header
module load R
echo "Running from $(pwd)"

R CMD BATCH --no-restore --no-save --quiet Flic2Sleep_script.R Flic_job.out