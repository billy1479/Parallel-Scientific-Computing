#!/bin/bash

# Request resources:
#SBATCH -c 15           # 1 CPU core
#SBATCH --mem=1G       # memory required, up to 250G on standard nodes.
#SBATCH --time=1:0:0   # time limit for job (format:  days-hours:minutes:seconds)
#SBATCH --gres=tmp:1G  # temporary disk space required on the compute node ($TMPDIR),
#                        up to 400G
# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared

# Commands to be run:
module load gcc
time ./10-2.sh