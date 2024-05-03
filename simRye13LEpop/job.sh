#!/bin/bash
#SBATCH -p ghpc_v3                # Name of the queue
#SBATCH -J hybrimax               # job name
#SBATCH --mem=8G                 # Memory in MiB(10 GiB = 10 * 1024 MiB)
#SBATCH -t 1000:00:00              # Job max time - Format = MM or MM:SS or HH:MM:SS or DD-HH or DD-HH:MM
#SBATCH -N 1                      # number of nodes
#SBATCH -n 1                      # number of cores
#SBATCH --output=slurm_%A.out     # STDOUT
#SBATCH --error=slurm_%A.err      # STDERR
#SBATCH --export=ALL
 
echo Job started at $(date '+%d_%m_%y_%H_%M_%S')
echo $SLURM_JOB_NODELIST
echo Job submitted from $SLURM_SUBMIT_HOST
echo Submission directory $SLURM_SUBMIT_DIR
 
ulimit -s unlimited
echo number of cores assigned: $SLURM_NTASKS
export OMP_NUM_THREADS=$SLURM_NTASKS

#find . -type f -name "slurm*.*" -delete
#find . -type f -name "out.log" -delete
#find . -type f -name "*.Rout" -delete

echo Running
R CMD BATCH --vanilla S4run1.R

#rdmuai -m 20480 -w 50:00:00 Inbredpop4Mod102

echo completed

#find . -type f -name "slurm*.*" -delete

#=========================================================================#
#                  Cleanup:  DO NOT REMOVE OR CHANGE                      #
#=========================================================================#
echo Job completed at $(date '+%d_%m_%y_%H_%M_%S')
