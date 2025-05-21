#!/bin/bash
# This line is required to inform the Linux
#command line to parse the script using
#the bash shell

# Instructing SLURM to locate and assign resources
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p "cpu"
#SBATCH -t 00-03:00:00

# Source the bash profile (required to use the module command)
source /etc/profile
source .venv/bin/activate

# Run your program (replace this with your program)
time ./step-5-instrument.sh