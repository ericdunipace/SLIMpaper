#!/bin/bash
#SBATCH -J Seed #A single job name for the array
#SBATCH -e error_seed.txt #Standard error
#SBATCH -p serial_requeue #Partition
#SBATCH -t 0-00:30:00 #Running time of 30 minutes
#SBATCH --mem-per-cpu 1000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH --mail-type=ALL #mail when start and finish
#SBATCH --mail-user=YOUR_EMAIL #set to your email for job notifications

module load gcc/7.1.0-fasrc01 R/3.5.1-fasrc02

source variables.txt

R CMD BATCH --no-save --no-restore "--args ${families[*]} cc ${corrs[*]} cc ${numobs[*]} cc ${numcoef[*]} cc ${penalty_factors[*]}" seed.R Output/seed.txt
