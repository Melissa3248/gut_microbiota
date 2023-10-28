#!/bin/bash
#SBATCH -N 1
#SBATCH -J sensitivity
#SBATCH -t 01:00:00
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

julia generate_params.jl 50 1 "/path/to/samples"