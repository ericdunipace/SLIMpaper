#!/bin/bash
#SBATCH -J Sparse-Posterior #A single job name for the array
#SBATCH -e SlurmOutput/sp_%A_%a_out.txt #Standard error
#SBATCH -o SlurmOutput/sp_%A_%a_out.txt #Standard output
#SBATCH -p serial_requeue #Partition
#SBATCH -t 0-16:30:00 #Running time of 16 hour(s) 30 min
#SBATCH --mem-per-cpu 5000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH --mail-type=ALL #mail when start and finish
#SBATCH --mail-user=YOUR_EMAIL #set to your email for job notifications

#module load gcc/7.1.0-fasrc01 R/3.5.1-fasrc02
module load gcc/8.2.0-fasrc01 openmpi/3.1.1-fasrc01 R/3.6.1-fasrc01

source variables.txt

R CMD BATCH --no-save "--args ${transport_methods[*]} cc ${penalty_factors[*]} cc ${PRIOR_SIGMA[*]} cc ${penalties[*]} cc ${families[*]} cc ${STAN_PATH[*]}" cluster.R TextOutput/${FAM}/${PENALTY}/${PENALTY_FACTOR}/${TRANSPORT_METHOD}/Corr_${CORR_X}/${NOBS}/${NCOEF}/SP_${TRANSPORT_METHOD}_Corr_${CORR_X}_${NOBS}_${NCOEF}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt
#R CMD BATCH --no-save "--args ${methods[*]} cc ${penalty_factors[*]} cc ${PRIOR_SIGMA[*]} cc ${penalties[*]} cc ${families[*]} cc ${x_methods[*]} cc ${STAN_PATH[*]}" cluster.R TextOutput/${X_METHOD}/${PENALTY}/${PENALTY_FACTOR}/${METHOD}/${NOBS}/${NCOEF}/SP_${METHOD}_${NOBS}_${NCOEF}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt
