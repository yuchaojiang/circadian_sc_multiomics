#!/bin/bash
#SBATCH --job-name=HMP       # Sets the job name
#SBATCH --time=14:00:00          # Sets the runtime
#SBATCH --ntasks=1              # Requests core
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G                # Requests memory
#SBATCH --out=HMP_%A-%a.out
#SBATCH --array=1-60

ml GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2
Rscript --no-save --no-restore --verbose 1_Resolution_tuning.R $SLURM_ARRAY_TASK_ID > HMP_$SLURM_ARRAY_TASK_ID.Rout 2>&1