'''
type this in the terminal

python run_parallel_solves.py  --save_directory "/path/to/save/directory"
'''

import os
import argparse

temp_job_folder = 'parallel_jobs'
def mkdir_p(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)
mkdir_p(temp_job_folder)

def run_job(temp_job_folder, job_name, num_samples, idx, save_directory):
    mkdir_p('{}'.format(temp_job_folder))
    job_file = '{}/{}.sh'.format(temp_job_folder, job_name)
    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash")
        fh.writelines("\n#SBATCH -N 1")
        fh.writelines("\n#SBATCH -C cpu")
        fh.writelines("\n#SBATCH -q regular # this is probably named something different on the iowa cluster")
        fh.writelines("\n#SBATCH -J {}".format(job_name))
        fh.writelines("\n#SBATCH -t 00:07:00")
        fh.writelines("\n#SBATCH --account=m4212")
        fh.writelines("\n#SBATCH --output={}/out_files/%x.out".format(temp_job_folder))
        fh.writelines("\njulia generate_params.jl {} {} {}".format(num_samples, idx, save_directory))
    os.system('sbatch {}'.format(job_file)) ################### uncomment this when you actually want to submit the jobs

######################################################################################
#                             Set up ArgumentParser                                  #
######################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("--save_directory", type=str)
args = parser.parse_args()

mkdir_p(args.save_directory)

ts = []

n_solves = 2048

n_runs = int(2**20/n_solves)

for t in range(n_runs):
    ts.append([n_solves,t, args.save_directory])


for t_ix in ts:
    num_samples = t_ix[0]
    idx = t_ix[1]
    save_directory = t_ix[2]
    job_name = "param_set{}".format(t_ix[1])
    run_job(temp_job_folder, job_name, num_samples, idx, save_directory)
