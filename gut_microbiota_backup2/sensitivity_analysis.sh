#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J sensitivity
#SBATCH --mail-user=maadrian@uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH -t 01:00:00
#SBATCH --account=m4212
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

julia generate_params.jl 50 1 "/pscratch/sd/m/maadrian/samples"