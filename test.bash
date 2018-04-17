#!/bin/bash
#SBATCH --time=200:00:00   
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p Long
. /etc/profile.d/modules.sh
module load matlab/R2016a-f

cd /home/xzh323/Metab
matlab -nodisplay -nodesktop -r 'run /home/xzh323/Metab/MCMCMetab/test.m'
wait
