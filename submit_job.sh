#!/bin/bash -login

### define resources needed:
### walltime - how long you expect the job to run
#SBATCH --time=5:30:00

### nodes:ppn - how many nodes & cores per node (ppn) that you require
#SBATCH --ntasks=2
### --constraint="lac"
#SBATCH --account="ptg"

### mem: amount of memory that the job will need
#SBATCH --mem-per-cpu=5G
### you can give your job a name for easier identification
#SBATCH -J pygmo_test 

### error/output file specifications
#SBATCH -o /mnt/home/herman67/cosy/pygmo/slurmfiles/moead_4f_2n.txt
#SBATCH -e /mnt/home/herman67/cosy/pygmo/slurmfiles/moead_4f_2n.txt
### load necessary modules, e.g.
#SBATCH --mail-user=herman67@msu.edu
#SBATCH --mail-type=FAIL
module restore gpflow
cd /mnt/home/herman67/cosy/pygmo
/mnt/home/herman67/anaconda3/envs/pygmo/bin/python optimize.py
scontrol show job ${SLURM_JOB_ID}

