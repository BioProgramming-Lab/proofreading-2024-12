#!/bin/bash

#SBATCH -p intel-sc3,amd-ep2,amd-ep2-short
#SBATCH -q huge
#SBATCH -J proofreading_simulation		
#SBATCH -c 1
#SBATCH -a 1-800
#SBATCH -o log/array-job.%A.%a.log 	

python figure_2_lhs.py