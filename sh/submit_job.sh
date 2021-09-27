#!/bin/bash -login

### define resources needed:
### walltime - how long you expect the job to run
#SBATCH --time=6-00:00:00
#SBATCH --begin=now

### nodes:ppn - how many nodes & cores per node (ppn) that you require
#SBATCH --ntasks=1
### --constraint="lac"
###SBATCH --account="ptg"

### mem: amount of memory that the job will need
#SBATCH --mem-per-cpu=5G
### you can give your job a name for easier identification
#SBATCH -J pygmo_test 
#SBATCH --array=120-129

### error/output file specifications
#SBATCH -o /mnt/simulations/secarml/secar_mobo/sh/slurmfiles/moead_4f_FP2_FP3_%a.txt
#SBATCH -e /mnt/simulations/secarml/secar_mobo/sh/slurmfiles/moead_4f_FP2_FP3_%a.txt
### load necessary modules, e.g.
###SBATCH --mail-user=herman67@msu.edu
###SBATCH --mail-type=FAIL
###module restore gpflow
cd /mnt/simulations/secarml/secar_mobo/py
./optimize.py ${SLURM_ARRAY_TASK_ID}
scontrol show job ${SLURM_JOB_ID}
