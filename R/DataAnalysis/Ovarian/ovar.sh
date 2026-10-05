#!/bin/bash
#SBATCH -J Ovarian_Selection #A single job name for the array
#SBATCH -e ../SlurmOutput/ovar_%A_%a_out.txt #Standard error
#SBATCH -o ../SlurmOutput/ovar_%A_%a_out.txt #Standard output
#SBATCH -p serial_requeue #Partition
#SBATCH -t 0-06:00:00 #Running time of 3 hour(s) 30 min
#SBATCH --mem-per-cpu 1000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH --mail-type=ALL #mail when start and finish
#SBATCH --mail-user=YOUR_EMAIL #set to your email for job notifications

#module load gcc/7.1.0-fasrc01 R/3.5.1-fasrc02
module load gcc/8.2.0-fasrc01 openmpi/3.1.1-fasrc01 R/3.6.1-fasrc01

R CMD BATCH --no-save ovarian_cluster_var.R ../TextOutput/Ovarian/Selection/ovarsel_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt
