#!/bin/bash
#SBATCH --array=1-25
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --job-name=minspan_call_100k

python ibd_sims/test/minspan_sweep.py call minspan_sweep_100k/constant_Ne_100k__coalescent \
    --iter $SLURM_ARRAY_TASK_ID