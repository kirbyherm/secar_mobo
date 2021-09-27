#!/bin/bash -login

### define resources needed:
### walltime - how long you expect the job to run
#SBATCH --time=0-03:30:00

### nodes:ppn - how many nodes & cores per node (ppn) that you require
#SBATCH --ntasks=1
### --constraint="lac"
###SBATCH --account="ptg"

### mem: amount of memory that the job will need
#SBATCH --mem-per-cpu=5G
### you can give your job a name for easier identification
#SBATCH -J pygmo_test 
###SBATCH --array=20-39

### error/output file specifications
#SBATCH -o /mnt/simulations/secarml/secar_mobo/sh/slurmfiles/moead_4f_FP2_FP3_draw.txt
#SBATCH -e /mnt/simulations/secarml/secar_mobo/sh/slurmfiles/moead_4f_FP2_FP3_draw.txt
### load necessary modules, e.g.
###SBATCH --mail-user=herman67@msu.edu
###SBATCH --mail-type=FAIL
###module restore gpflow
cd /mnt/simulations/secarml/secar_mobo/py/4d_problem/4f_FP2_FP3
for TOSSES in $(seq 28 1 28) 
do
    cosy pygmoCosy${TOSSES}.fox
    mv pic001.pdf X${TOSSES}.pdf
    mv pic002.pdf Y${TOSSES}.pdf
done

#scontrol show job ${SLURM_JOB_ID}

