'''
type this in the terminal

python run_parallel_solves_UI.py  --save_directory "/path/to/save/directory"

#!/bin/bash

#####Set Scheduler Configuration Directives#####
#Set the name of the job. This will be the first part of the error/output filename.
#$ -N {}

#Set the current working directory as the location for the error and output files.
#(Will show up as .e and .o files)
#$ -cwd {}

#####End Set Scheduler Configuration Directives#####


#####Resource Selection Directives#####
#See the HPC wiki for complete resource information: https://wiki.uiowa.edu/display/hpcdocs/Argon+Cluster
#Select the queue to run in
#$ -q BA

#Select the number of slots the job will use # number of threads for the process
#$ -pe smp 1 

#####Begin Compute Work#####
#Print information from the job into the output file
/bin/echo Running on compute node: `hostname`.
/bin/echo In directory: `pwd`
/bin/echo Starting on: `date`

#Sleep (wait and do nothing) for 60 seconds
#module load stack/2020.1
module load julia
julia generate_params.jl {} {} {}

#Print the end date of the job before exiting
echo Now it is: `date`
#####End Compute Work#####



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
    job_file = '{}/{}.job'.format(temp_job_folder, job_name)
    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash")
        fh.writelines("\n#$ -N {}".format(job_name)) # name of the job
        fh.writelines("\n#$ -q BA") # name of the queue
        fh.writelines("\n#$ -pe smp 1") # number of threads to use
        #fh.writelines("\n#SBATCH -q regular # this is probably named something different on the iowa cluster")
        #fh.writelines("\n#SBATCH -J {}".format(job_name))
        #fh.writelines("\n#SBATCH -t 00:15:00")
        #fh.writelines("\n#SBATCH --account=m4212_g")
        fh.writelines("\n#$ -cwd") # Determines whether the job will be executed from the current working directory. 
                                    # If not specified, the job will be run from your home directory.
        fh.writelines("\n#$ -e {}/out_files/{}.out".format(temp_job_folder,job_name)) #) # directory to store the out files
        fh.writelines("\n#module load stack/2020.1")
        fh.writelines("module load julia")
        fh.writelines("\njulia generate_params.jl {} {} {}".format(num_samples, idx, save_directory))
    #os.system('qsub {}'.format(job_file)) ################### uncomment this when you actually want to submit the jobs

######################################################################################
#                             Set up ArgumentParser                                  #
######################################################################################
parser = argparse.ArgumentParser()
parser.add_argument("--save_directory", type=str)
args = parser.parse_args()

mkdir_p(args.save_directory)

ts = []

n_runs = int(2**20/512)

for t in range(n_runs):
    ts.append([512,t, args.save_directory])


for t_ix in ts:
    num_samples = t_ix[0]
    idx = t_ix[1]
    save_directory = t_ix[2]
    job_name = "param_set{}".format(t_ix[1])
    run_job(temp_job_folder, job_name, num_samples, idx, save_directory)
