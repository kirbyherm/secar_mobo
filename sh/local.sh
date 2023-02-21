#!/bin/bash
SLURM_ARRAY_TASK_ID=1

cd ~/code/secar_mobo/py
./optimize.py ${SLURM_ARRAY_TASK_ID}

