#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular # this is probably named something different on the iowa cluster
#SBATCH -J param_set993
#SBATCH -t 00:15:00
#SBATCH --output=parallel_jobs/out_files/%x.out
julia generate_params.jl 512 993 /pscratch/sd/m/maadrian/testing_saving