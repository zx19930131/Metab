#!/bin/bash
#SBATCH --time=200:00:00   
#SBATCH -N 1
#SBATCH -n 1
cd /home/xzh323/Metab/MCMCMetab
R CMD BATCH test.m
wait  
