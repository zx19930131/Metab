#!/bin/bash
#SBATCH --time=200:00:00   
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p Long
cd /home/xzh323/Metab
module load matlab/R2016a-f
matlab -nodisplay -nodesktop -r 'run /home/xzh323/Metab/MCMCMetab/test4.m'
wait  
