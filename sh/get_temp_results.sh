#!/bin/bash -login 
### define resources needed:
### walltime - how long you expect the job to run
#SBATCH --time=0-02:11:00
#SBATCH --begin=now

### nodes:ppn - how many nodes & cores per node (ppn) that you require
###SBATCH --ntasks=10
###SBATCH --nodes=1
### --constraint="lac"
###SBATCH --account="ptg"

### mem: amount of memory that the job will need
#SBATCH --mem-per-cpu=5G
### you can give your job a name for easier identification
#SBATCH -J cpu_test 
#SBATCH --array=19

### error/output file specifications
#SBATCH -o /mnt/simulations/secarml/secar_mobo/sh/slurmfiles/test_%a.txt
#SBATCH -e /mnt/simulations/secarml/secar_mobo/sh/slurmfiles/test_%a.txt
### load necessary modules, e.g.
###SBATCH --mail-user=herman67@msu.edu
###SBATCH --mail-type=FAIL
###module restore gpflow
#cd /mnt/simulations/secarml/secar_mobo/py
#./cosy.py
echo "Checking scratch"
cd /scratch/
echo "Checking hermanse"
cd hermanse
ls -lrt
echo "Checking hermanse/secar_mobo/"
cd secar_mobo
ls -lrt
#rm -rf 10${SLURM_ARRAY_TASK_ID}
cd 20${SLURM_ARRAY_TASK_ID}
ls -lrt
cat output20${SLURM_ARRAY_TASK_ID}.csv > /mnt/simulations/secarml/secar_mobo/output/output_temp20${SLURM_ARRAY_TASK_ID}.csv
scontrol show job ${SLURM_JOB_ID}
#hostname
#hostnamectl
