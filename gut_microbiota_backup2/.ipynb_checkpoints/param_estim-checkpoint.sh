#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J param_estim
#SBATCH --mail-user=maadrian@uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH -t 01:00:00
#SBATCH --account=m4212

python param_estim.py