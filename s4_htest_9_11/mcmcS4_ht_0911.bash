#!/bin/bash
#SBATCH --time 14-00:00:00     # Time limit for the job (REQUIRED)
#SBATCH --job-name=ht0911    # Job name
#SBATCH --nodes=1        # Number of nodes to allocate
#SBATCH --ntasks=1       # Number of cores to allocate
#SBATCH --partition=SKY32M192_L     # Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err  # Error file
#SBATCH -o slurm-%j.out  # Output file
#SBATCH -A col_cwa236_uksr  # Project allocation account name (REQUIRED)
#SBATCH --mail-type ALL    # Send email when job starts/ends
#SBATCH --mail-user xzh323@g.uky.edu   # Where email is sent to (optional)

module load gnu/5.4.0
module ccs/boost/1.70.0
module ccs/eigen/3.3.5

cd /home/xzh323/Metab/s4_htest_9_11
./myMCMC
wait
