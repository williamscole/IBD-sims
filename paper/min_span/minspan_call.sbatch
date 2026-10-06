#!/bin/bash
#SBATCH --job-name=minspan_call
#SBATCH --array=1-25
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=minspan_sweep/logs/minspan_call_%A_%a.log

# activate your environment here
python ibd_sims/test/minspan_sweep.py call minspan_sweep/constant_Ne_10k__coalescent --iter $SLURM_ARRAY_TASK_ID
