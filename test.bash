#!/bin/bash
#SBATCH --time=200:00:00   
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p Long
cd /home/xzh323/Metab/MCMCMetab
matlab -nodisplay -nodesktop -r /home/xzh323/Metab/MCMCMetab/test.m
wait  
