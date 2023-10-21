#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular # this is probably named something different on the iowa cluster
#SBATCH -J param_set1883
#SBATCH -t 00:15:00
#SBATCH --output=parallel_jobs/out_files/%x.out
julia generate_params.jl 512 1883 /pscratch/sd/m/maadrian/testing_saving